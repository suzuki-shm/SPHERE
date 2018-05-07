#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-11-01

from sphere.sphere_utils import compress_depth
from sphere.sphere_utils import load_depth_file
from sphere.sphere_utils import get_logger
import argparse
import pandas as pd
import numpy as np


def argument_parse(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("depth_file_path",
                        type=str,
                        help="file path of coverage depth")
    parser.add_argument("output_dest",
                        type=str,
                        help="destination of output tsv file")
    parser.add_argument("-cl", "--compressedlength",
                        dest="cl",
                        nargs="?",
                        default=10000,
                        type=int,
                        help="Compressed length of genome (default: 10000)")
    args = parser.parse_args(argv)
    return vars(args)


def main(args, logger):
    df = load_depth_file(args["depth_file_path"])

    genome_name = df["genome"].unique()[0]
    position = np.arange(1, args["cl"]+1, 1)
    c_depth = compress_depth(df["depth"].values, args["cl"])
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
