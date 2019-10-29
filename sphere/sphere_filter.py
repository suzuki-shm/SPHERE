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
                            "sum",
                            "comp"
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
                        help="Allowed range of variance filter (default: 2)")
    parser.add_argument("-p", "--percentile",
                        dest="p",
                        nargs="?",
                        default=0.99,
                        type=float,
                        help="Threshold of percentile filter (default: 0.99)")
    parser.add_argument("-l", "--length",
                        dest="l",
                        nargs="?",
                        default=None,
                        type=int,
                        help="Length to be compensate (default: None)")
    parser.add_argument("-m", "--missing",
                        dest="m",
                        nargs="?",
                        default=np.nan,
                        help="Filling value by the filter (default: nan)")
    args = parser.parse_args(argv)
    return vars(args)


def main(args, logger):
    df = load_depth_file(args["depth_file_path"])
    if args["t"] in ["median", "sum", "variance"]:
        # Evaluate input parameter
        if len(df) < args["w"]:
            raise ValueError("Window size must be smaller than input size")
        elif len(df) < args["s"]:
            raise ValueError("Stride size must be smaller than input size")

        f_df = df.groupby(["genome"])["depth"].apply(moving_filter,
                                                     s=args["s"],
                                                     w=args["w"],
                                                     r=args["r"],
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
    elif args["t"] == "percentile":
        percentile_value = df.groupby("genome")["depth"].quantile(args["p"])
        index = df[df["depth"] > df["genome"].apply(lambda x: percentile_value[x])].index
        index_nzero = df[df["depth"] > 0].index
        if len(index) == len(df):
            msg = " ".join([
                "The filtering process aborted ",
                "as this percentile will remove all coverage depth"
            ])
            logger.warning(msg)
        elif len(index_nzero) == len(index):
            msg = " ".join([
                "The filterin process aborted ",
                "as this percentile will convert all depth to zeros"
            ])
            logger.warning(msg)
        else:
            df.loc[index, "depth"] = args["m"]
        f_df = df
    elif args["t"] == "fill":
        f_df = df
        if args["m"] == "median":
            f_df["depth"] = f_df["depth"].fillna(f_df["depth"].median())
        elif args["m"] == "nan":
            f_df["depth"] = f_df["depth"].fillna(np.nan)
        elif args["m"] == "mmedian":
            f_df_mm = f_df["depth"].rolling(window=args["w"],
                                            min_periods=1,
                                            center=True).median()
            f_df_nan_index = f_df[f_df["depth"].isna()].index
            f_df.loc[f_df_nan_index, "depth"] = f_df_mm[f_df_nan_index]
            if f_df["depth"].isna().sum() != 0:
                raise ValueError("The result still contains NaN. "
                                 "Use more greater window size.")
        else:
            f_df["depth"] = f_df["depth"].fillna(args["m"])
        f_df["depth"] = f_df["depth"].astype(int)
    elif args["t"] == "comp":
        if len(df["genome"].unique()) != 1:
            raise ValueError(
                "For comp mode,",
                "the input file must be aligned to single reference sequence."
            )
        genome = df["genome"].unique()[0]
        f_df = pd.DataFrame(
            {"location": range(1, args["l"]+1)}
        )
        f_df = f_df.merge(df, on="location", how="left")
        f_df["genome"] = genome
        f_df["depth"] = f_df["depth"].fillna(0)
        f_df["depth"] = f_df["depth"].astype(int)
        f_df = f_df[["genome", "location", "depth"]]

    f_df.to_csv(args["output_dest"], sep="\t", index=None, header=None)


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
