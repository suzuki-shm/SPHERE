#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-11-01

from sphere.sphere_utils import moving_filter
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
                        choices=[
                            "median",
                            "variance",
                            "percentile",
                            "fill",
                            "sum"
                        ],
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
                        help="Filling value by the filter (default: nan)")
    args = parser.parse_args(argv)
    return vars(args)


def main(args, logger):
    df = load_depth_file(args["depth_file_path"])
    if args["t"] == "median" or args["t"] == "sum":
        f_df = df.groupby(["genome"])["depth"].apply(moving_filter,
                                                     s=args["s"],
                                                     w=args["w"],
                                                     ftype=args["t"])
        # Need if statement
        # because of bug in pandas
        # (https://github.com/pandas-dev/pandas/issues/13056)
        if type(f_df.index) is pd.MultiIndex:
            f_df.index.set_names("location", level=1, inplace=True)
        else:
            f_df.index.set_names("location", inplace=True)
        f_df = f_df.reset_index()
        f_df["location"] += 1
        if f_df.shape[1] == 2:
            f_df["genome"] = df["genome"][0]
        f_df = f_df[["genome", "location", "depth"]]
    elif args["t"] == "variance":
        p = 1 / df.groupby("genome")["location"].max()
        n = df.groupby("genome")["depth"].sum()
        m = df.groupby("genome")["depth"].mean()
        v = np.sqrt(n * p * (1-p))
        m_dict = m.to_dict()
        v_dict = v.to_dict()
        depth_diff = abs(df["genome"].apply(lambda x: m_dict[x]) - df["depth"])
        depth_std = args["r"] * df["genome"].apply(lambda x: v_dict[x])
        index = df[depth_diff < depth_std].index
        df.loc[index, "depth"] = args["m"]
        f_df = df
    elif args["t"] == "percentile":
        percentile_value = df.groupby("genome")["depth"].quantile(args["p"])
        index = df[df["depth"] > df["genome"].apply(lambda x: percentile_value[x])].index
        df.loc[index, "depth"] = args["m"]
        f_df = df
    elif args["t"] == "fill":
        f_df = df
        if args["m"] == "median":
            f_df["depth"] = f_df["depth"].fillna(f_df["depth"].median())
        elif args["m"] == "nan":
            f_df["depth"] = f_df["depth"].fillna(np.nan)
        else:
            f_df["depth"] = f_df["depth"].fillna(args["m"])
        f_df["depth"] = f_df["depth"].astype(int)

    f_df.to_csv(args["output_dest"], sep="\t", index=None, header=None)


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
