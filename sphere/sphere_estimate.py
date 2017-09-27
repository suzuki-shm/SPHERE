#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

from logging import getLogger, DEBUG, Formatter, StreamHandler
from sphere.sphere_compile import compile_model
from sphere.stan_utils import save_model
from sphere.stan_utils import load_model
from sphere.stan_utils import summarize_fit
from sphere.stan_utils import save_fit
from sphere.stan_utils import sampling
from sphere.sphere_utils import compress_depth
from sphere.sphere_utils import load_depth_file
import os
import argparse


def get_logger():
    logger = getLogger(__name__)
    logger.setLevel(DEBUG)
    log_fmt = '%(asctime)s : %(name)s : %(levelname)s : %(message)s'
    formatter = Formatter(log_fmt)
    stream_handler = StreamHandler()
    stream_handler.setLevel(DEBUG)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    return logger


def argument_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("depth_file_path",
                        type=str,
                        help="file path of coverage depth")
    parser.add_argument("output_dest",
                        type=str,
                        help="destination of output tsv file")
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
    parser.add_argument("-cl", "--compressedlength",
                        dest="cl",
                        nargs="?",
                        default=10000,
                        type=int,
                        help="Compressed length of genome (default: 10000)")
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
    parser.set_defaults(trans=False)
    args = parser.parse_args()
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
    df = load_depth_file(args["depth_file_path"])
    if args["pmp"] is not None:
        model = load_model(args["compiled_model_path"])
    else:
        model = compile_model(args["pmd"])
    if args["pmd"] is not None:
        save_model(args["pmd"], model)
    I = len(df)
    v_c = compress_depth(df["depth"], I, args["cl"])
    fit = sampling(model,
                   v_c,
                   args["si"], args["sw"], args["sc"], args["st"], args["ss"])
    sdf = summarize_fit(fit)
    sdf.to_csv(args["output_dest"], sep="\t")
    if args["fod"] is not None:
        save_fit(fit, args["fod"])


def main_wrapper():
    args = argument_parse()
    logger = get_logger()
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
