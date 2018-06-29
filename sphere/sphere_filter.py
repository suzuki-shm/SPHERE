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
                        help="destination of output tsv file")
    parser.add_argument("depth_file_path",
                        type=str,
                        help="file path of coverage depth")
    parser.add_argument("-s", "--stride_length",
                        dest="s",
                        nargs="?",
                        default=100,
                        type=int,
                        help="Stride length of filter (default: 100)")
    parser.add_argument("-w", "--window_length",
                        dest="w",
                        nargs="?",
                        default=10000,
                        type=int,
                        help="Window length of filter (default: 10000)")
    args = parser.parse_args(argv)
    return vars(args)


def main(args, logger):
    df = load_depth_file(args["depth_file_path"])
    cl = compress_length(df["depth"].size, s=args["s"], w=args["w"])

    genome_name = df["genome"].unique()[0]
    position = np.arange(1, cl + 1, 1).astype(int)
    c_depth = compress_depth(df["depth"], s=args["s"], w=args["w"])
    c_df = pd.DataFrame({"position": position, "depth": c_depth})
    c_df["genome"] = genome_name
    c_df = c_df[["genome", "position", "depth"]]

    c_df.to_csv(args["output_dest"], sep="\t", index=None, header=None)


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
