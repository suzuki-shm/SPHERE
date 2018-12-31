#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-11-01

from sphere.sphere_utils import compress_depth
from sphere.sphere_utils import compress_length
from sphere.sphere_utils import load_depth_file
from sphere.sphere_utils import get_logger
import argparse
import pandas as pd
import numpy as np


def argument_parse(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dest",
                        type=str,
                        help="Destination of output tsv file")
    parser.add_argument("depth_file_path",
                        type=str,
                        help="File path of coverage depth")
    parser.add_argument("-t", "--type",
                        dest="t",
                        nargs="?",
                        default="median",
                        choices=["median", "variance", "percentile", "fill"],
                        type=str,
                        help="Filter type (default: median)")
    parser.add_argument("-s", "--stride_length",
                        dest="s",
                        nargs="?",
                        default=100,
                        type=int,
                        help="Stride length of median filter (default: 100)")
    parser.add_argument("-w", "--window_length",
                        dest="w",
                        nargs="?",
                        default=10000,
                        type=int,
                        help="Window length of median filter (default: 10000)")
    parser.add_argument("-r", "--range",
                        dest="r",
                        nargs="?",
                        default=2,
                        type=int,
                        help="Range of variance filter (default: 2)")
    parser.add_argument("-p", "--percentile",
                        dest="p",
                        nargs="?",
                        default=0.99,
                        type=float,
                        help="Threshold of percentile filter (default: 0.99)")
    parser.add_argument("-m", "--missing",
                        dest="m",
                        nargs="?",
                        default=np.nan,
                        help="Filling value by the filter (default: 99999)")
    args = parser.parse_args(argv)
    return vars(args)


def main(args, logger):
    df = load_depth_file(args["depth_file_path"])
    if args["t"] == "median":
        cl = compress_length(df["depth"].size, s=args["s"], w=args["w"])
        genome_name = df["genome"].unique()[0]
        location = np.arange(1, cl + 1, 1).astype(int)
        c_depth = compress_depth(df["depth"], s=args["s"], w=args["w"])
        f_df = pd.DataFrame({"location": location, "depth": c_depth})
        f_df["genome"] = genome_name
        f_df = f_df[["genome", "location", "depth"]]
    elif args["t"] == "variance":
        p = df["location"].max()
        n = df["depth"].sum()
        m = df["depth"].mean()
        v = n * p * (1-p)
        index = df.query("depth - {0} > {0}*{1}".format(m, args["r"], v)).index
        df.loc[index, "depth"] = args["m"]
        f_df = df
    elif args["t"] == "percentile":
        percentile_value = df["depth"].quantile(args["p"])
        index = df.query("depth > {0}".format(percentile_value)).index
        df.loc[index, "depth"] = args["m"]
        f_df = df
    elif args["t"] == "fill":
        f_df = df.fillna(args["m"])
        f_df["depth"] = f_df["depth"].astype(int)

    f_df.to_csv(args["output_dest"], sep="\t", index=None, header=None)


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
