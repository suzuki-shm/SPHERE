#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

from sphere.sphere_compile import compile_model
from sphere.stan_utils import save_model
from sphere.stan_utils import load_model
from sphere.stan_utils import summarize_fit
from sphere.stan_utils import save_fit
from sphere.stan_utils import sampling
from sphere.sphere_utils import get_logger
from sphere.sphere_utils import load_multiple_depth_file
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
    parser.add_argument("-pmd", "--pickled_model_dest",
                        dest="pmd",
                        nargs="?",
                        default=None,
                        type=str,
                        help="destination of pickled model (default: None)")
    parser.add_argument("-pmp", "--pickled_model_path",
                        dest="pmp",
                        nargs="?",
                        default=None,
                        type=str,
                        help="file path of compiled model (default: None)")
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
                        default="trigonal",
                        type=str,
                        choices=["trigonal",
                                 "linear",
                                 "vonmises"],
                        help="model type for trend (default: trigonal)")
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
                        default=3,
                        type=int,
                        help="Number of Stan chain (default: 3)")
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
    parser.set_defaults(ff=False)
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

    if args["pmp"] is not None:
        model = load_model(args["pmp"])
    else:
        logger.info("Compiling stan model")
        model = compile_model(args["pmd"], args["m"])
    if args["pmd"] is not None:
        logger.info("Saving compiled model to {0}".format(args["pmd"]))
        save_model(args["pmd"], model)

    stan_data = {
        "I": len(df),
        "S": len(args["depth_file_path"]),
        "L": df["location"].max(),
        "SUBJECT": df["subject"].values,
        "LOCATION": df["location"].values,
        "DEPTH": df["depth"].values
    }
    logger.info("Sampling from probability distribution")
    fit = sampling(model,
                   stan_data,
                   args["p"],
                   args["si"], args["sw"], args["sc"], args["st"], args["ss"])
    logger.info("Summarizing MCMC result")
    sdf = summarize_fit(fit, pars=args["p"])
    logger.info("Saving MCMC summary to {0}".format(args["output_dest"]))
    sdf.to_csv(args["output_dest"], sep="\t")
    if args["fod"] is not None:
        logger.info("Saving MCMC all result to {0}".format(args["fod"]))
        save_fit(fit, args["fod"])
    if args["lld"] is not None:
        logger.info("Saving log likelifood to {0}".format(args["lld"]))
        save_log_lik(fit, args["lld"])


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
