#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

from sphere.sphere_utils import get_logger
from sphere.sphere_utils import load_multiple_depth_file
import argparse
from pystan.misc import rdump


def argument_parse(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dest",
                        type=str,
                        help="destination of output tsv file")
    parser.add_argument("depth_file_path",
                        type=str,
                        nargs="+",
                        help="file pathes of coverage depth")
    parser.add_argument("-nmix", "--number_mixture",
                        dest="nmix",
                        nargs="?",
                        default=1,
                        type=int,
                        help="Number of mixed distribution (default: 1)")
    args = parser.parse_args(argv)
    return vars(args)


def main(args, logger):
    logger.info("Loading sequence depth file")
    df = load_multiple_depth_file(args["depth_file_path"])
    n_samples = len(args["depth_file_path"])
    n_length = df["location"].max()
    # Drop tuples those depth is 0 to reduce memory usage

    stan_data = {}
    stan_data["K"] = args["nmix"]
    df = df[df["depth"] != 0]
    n_iteration = len(df)
    stan_data["I"] = n_iteration
    stan_data["S"] = n_samples
    stan_data["L"] = n_length
    stan_data["SUBJECT"] = df["subject"].values
    stan_data["LOCATION"] = df["location"].values
    stan_data["DEPTH"] = df["depth"].values

    rdump(stan_data, args["output_dest"])


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
