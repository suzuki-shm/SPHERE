#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

from sphere.stan_utils import load_model
from sphere.stan_utils import summarize_fit
from sphere.stan_utils import summarize_ofit
from sphere.stan_utils import save_fit
from sphere.stan_utils import sampling
from sphere.stan_utils import optimizing
from sphere.sphere_utils import get_logger
from sphere.sphere_utils import load_multiple_depth_file
from sphere.sphere_utils import get_pars
from sphere.stan_utils import save_log_lik
import os
import argparse


def argument_parse(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dest",
                        type=str,
                        help="destination of output tsv file")
    parser.add_argument("depth_file_path",
                        type=str,
                        nargs="+",
                        help="file pathes of coverage depth")
    parser.add_argument("-fod", "--fit_out_dest",
                        dest="fod",
                        nargs="?",
                        default=None,
                        type=str,
                        help="file path of estimated result (default: None)")
    parser.add_argument("-lld", "--log_likelihood_dest",
                        dest="lld",
                        nargs="?",
                        default=None,
                        type=str,
                        help="file path of estimated lig_lik (default: None)")
    parser.add_argument("-m", "--model",
                        dest="m",
                        nargs="?",
                        default="vonmises",
                        type=str,
                        choices=[
                            "linearcardioid",
                            "cardioid",
                            "wrappedcauchy",
                            "vonmises",
                            "jonespewsey",
                            "dvonmises",
                            "miaecardioid",
                            "miaewrappedcauchy",
                            "miaevonmises",
                            "miaejonespewsey",
                            "invmiaecardioid",
                            "invmiaewrappedcauchy",
                            "invmiaevonmises",
                            "invmiaejonespewsey"
                        ],
                        help="model type for trend (default: vonmises)")
    parser.add_argument("-M", "--method",
                        dest="M",
                        nargs="?",
                        default="sampling",
                        type=str,
                        choices=["sampling",
                                 "optimizing"],
                        help="adaptation method")
    parser.add_argument("-nmix", "--number_mixture",
                        dest="nmix",
                        nargs="?",
                        default=1,
                        type=int,
                        help="Number of mixed distribution (default: 1)")
    parser.add_argument("-si", "--staniter",
                        dest="si",
                        nargs="?",
                        default=3000,
                        type=int,
                        help="Number of Stan iteration (default: 3000)")
    parser.add_argument("-sw", "--stanwarmup",
                        dest="sw",
                        nargs="?",
                        default=1000,
                        type=int,
                        help="Number of Stan warm up (defaultL 1000)")
    parser.add_argument("-sc", "--stanchain",
                        dest="sc",
                        nargs="?",
                        default=1,
                        type=int,
                        help="Number of Stan chain (default: 1)")
    parser.add_argument("-st", "--stanthin",
                        dest="st",
                        nargs="?",
                        default=1,
                        type=int,
                        help="Number of Stan thin (default: 1)")
    parser.add_argument("-ss", "--stanseed",
                        dest="ss",
                        nargs="?",
                        default=1234,
                        type=int,
                        help="Number of Stan seed (defaultt: 1234)")
    parser.add_argument("-ff",
                        dest="ff",
                        action="store_true",
                        help="overwrite output file (default: False)")
    parser.add_argument("-p", "--pars",
                        dest="p",
                        nargs="*",
                        default=None,
                        type=str,
                        help="parameter of interest (default: None)")
    parser.add_argument("-ll", "--log_lik",
                        dest="ll",
                        action="store_true",
                        help="log_lik contained in output?(default: False)")
    parser.add_argument("-j", "--n_jobs",
                        dest="j",
                        default=-1,
                        type=int,
                        help="Number of parallel chain used only in sampling")
    parser.set_defaults(ff=False, ll=False)
    args = parser.parse_args(argv)
    return vars(args)


def check_output_dest(ff, od):
    if ff is False:
        if os.path.exists(od):
            raise ValueError("{0} has been already existed. ".format(od),
                             "Use -ff to overwrite.")
        else:
            return 0
    else:
        if os.path.exists(od):
            return 1
        else:
            return 0


def main(args, logger):
    od_exist = check_output_dest(args["ff"], args["output_dest"])
    if od_exist:
        logger.warning("{0} will be overwrited".format(args["output_dest"]))
    logger.info("Loading sequence depth file")
    df = load_multiple_depth_file(args["depth_file_path"])
    n_samples = len(args["depth_file_path"])
    n_length = df["location"].max()
    # Drop tuples those depth is 0 to reduce memory usage

    stan_data = {}
    if args["sc"] != 1 and args["nmix"] > 1:
        msg = "As number of chains must be one for mixture model "
        msg += "to aboid label switching, it is setted to one."
        logger.warning(msg)
        args["si"] = args["si"] * args["sc"]
        args["sc"] = 1
    stan_data["K"] = args["nmix"]
    df = df[df["depth"] != 0]
    n_iteration = len(df)
    stan_data["I"] = n_iteration
    stan_data["S"] = n_samples
    stan_data["L"] = n_length
    stan_data["SUBJECT"] = df["subject"].values
    stan_data["LOCATION"] = df["location"].values
    stan_data["DEPTH"] = df["depth"].values

    logger.info("Loading model file")
    model = load_model(args["m"])
    if args["p"] is None:
        if args["lld"] is not None and args["ll"] is False:
            msg = "Outputting log likelihood of posterior distribution fails,"
            msg += "if ll parameter is not used."
            msg += "Adding log_lik parameter in MCMC sampling."
            logger.warning(msg)
            args["ll"] = True
        pars = get_pars(args["m"], args["ll"])
    else:
        pars = args["p"]

    if args["M"] == "sampling":
        logger.info("Sampling from probability distribution")
        fit = sampling(model,
                       stan_data,
                       pars,
                       args["si"], args["sw"], args["sc"], args["st"],
                       args["ss"], args["j"])
        logger.info("Summarizing result")
        sdf = summarize_fit(fit, pars=pars)
        logger.info("Saving summary to {0}".format(args["output_dest"]))
        sdf.to_csv(args["output_dest"], sep="\t")
        if args["fod"] is not None:
            logger.info("Saving MCMC all result to {0}".format(args["fod"]))
            save_fit(fit, args["fod"])
        if args["lld"] is not None:
            logger.info("Saving log likelifood to {0}".format(args["lld"]))
            save_log_lik(fit, args["lld"])
    elif args["M"] == "optimizing":
        logger.info("Optimizing the parameters to the data")
        try:
            ofit = optimizing(model, stan_data, args["ss"])
        except RuntimeError:
            logger.warning("L-BFGS algorithm failed. Use Newton algorithm")
            ofit = optimizing(model, stan_data, args["ss"], "Newton")
        logger.info("Summarizing result")
        sdf = summarize_ofit(ofit, pars=pars)
        logger.info("Saving summary to {0}".format(args["output_dest"]))
        sdf.to_csv(args["output_dest"], sep="\t")
    else:
        raise ValueError("Invalid argument for method")


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
