#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-28

from logging import getLogger, DEBUG, Formatter, StreamHandler
from sphere.sphere_utils import load_depth_file
import argparse
import numpy as np
import pandas as pd
import os


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
                        nargs="*",
                        type=str,
                        help="file(s) path of coverage depth")
    parser.add_argument("output_dest",
                        type=str,
                        help="destination of output tsv file")
    args = parser.parse_args()
    return vars(args)


def mean_resultant_length(C, S, y_sum):
    R = np.sqrt(C**2 + S**2) / y_sum
    return R


def circular_variance(R):
    V = 1 - R
    return V


def circular_standard_deviation(R):
    v = np.sqrt(-2 * np.log(R))
    return v


def mean_direction(S, C):
    if C > 0 and S >= 0:
        return np.arctan(S/C)
    elif C == 0 and S > 0:
        return np.pi/2
    elif C < 0:
        return np.arctan(S/C) + np.pi
    elif C == 0 and S < 0:
        return 3 * np.pi / 2
    else:
        return np.arctan(S/C) + 2 * np.pi


def mean_direction_1dimension(md, I):
    md_1d = md / (2 * np.pi) * I
    if md_1d < 0:
        return md_1d + I
    else:
        return md_1d


def main(args, logger):
    result = []
    for f in args["depth_file_path"]:
        file_name = os.path.basename(f)
        df = load_depth_file(f)
        I = len(df)
        x = np.arange(1, I+1, 1)
        y = df["depth"].values
        x_cos = np.cos(x / I * 2 * np.pi)
        x_sin = np.sin(x / I * 2 * np.pi)
        C = np.sum(y * x_cos)
        S = np.sum(y * x_sin)

        mrl = mean_resultant_length(C, S, np.sum(y))
        cv = circular_variance(mrl)
        csd = circular_standard_deviation(mrl)
        md = mean_direction(S, C)
        md_1d = mean_direction_1dimension(md, I)

        tmp = {
            "filepath": f,
            "filename": file_name,
            "mean_resultant_length": mrl,
            "circular_variance": cv,
            "circular_standard_deviation": csd,
            "mean_direction": md,
            "mean_position": md_1d
        }
        result.append(tmp)
    result_df = pd.DataFrame(result)
    result_df = result_df.set_index("filename")
    result_df = result_df[[
                            "mean_resultant_length",
                            "circular_variance",
                            "circular_standard_deviation",
                            "mean_direction",
                            "mean_position",
                            "filepath"
                         ]]
    result_df.to_csv(args["output_dest"], sep="\t")


def main_wrapper():
    args = argument_parse()
    logger = get_logger()
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
