#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2018-04-26

from sphere.sphere_utils import load_depth_file
from scipy.stats import norm
import argparse
import numpy as np
import pandas as pd
import os


def argument_parse(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dest",
                        type=str,
                        help="destination of output tsv file")
    parser.add_argument("depth_file_path",
                        nargs="+",
                        type=str,
                        help="file(s) path of coverage depth")
    args = parser.parse_args(argv)
    return vars(args)


def mean_resultant_length(C, S, y_sum):
    # y_sum means number of sample directions
    R = np.sqrt(C**2 + S**2) / y_sum
    return R


def mean_direction(S, C):
    if C > 0 and S >= 0:
        return np.arctan(S / C)
    elif C == 0 and S > 0:
        return np.pi / 2.0
    elif C < 0:
        return np.arctan(S/C) + np.pi
    elif C == 0 and S < 0:
        return 3.0 * np.pi / 2.0
    else:
        return np.arctan(S / C) + 2 * np.pi


def sin_moment(theta, y, md, p=1):
    return (y * np.sin(p * (theta - md))).mean()


def cos_moment(theta, y, md, p=1):
    return (y * np.cos(p * (theta - md))).mean()


def perwey_test(theta, y):
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    C = np.sum(y * cos_theta)
    S = np.sum(y * sin_theta)

    mrl = mean_resultant_length(C, S, np.sum(y))
    md = mean_direction(S, C)
    b2 = sin_moment(theta, y, md, p=2)
    a2 = cos_moment(theta, y, md, p=2)
    a3 = cos_moment(theta, y, md, p=3)
    a4 = cos_moment(theta, y, md, p=4)
    var_b2 = ((1 - a4)/2 - 2*a2 + 2*a2/mrl*(a3+(a2*(1-a2))/mrl))

    z = b2 / np.sqrt(var_b2)
    p = 1 - norm.cdf(abs(z))

    return (z, p)


def main(args):
    result = []
    for f in args["depth_file_path"]:
        file_name = os.path.basename(f)
        df = load_depth_file(f)
        I = len(df)
        x = np.arange(1, I+1, 1)
        theta = x / float(I) * 2.0 * np.pi
        y = df["depth"].values
        z, p = perwey_test(theta, y)

        tmp = {
            "filepath": f,
            "filename": file_name,
            "z": z,
            "p": p
        }
        result.append(tmp)
    result_df = pd.DataFrame(result)
    result_df = result_df.set_index("filename")
    result_df = result_df[[
                            "z",
                            "p",
                            "filepath"
                         ]]
    result_df.to_csv(args["output_dest"], sep="\t")


def main_wrapper():
    args = argument_parse()
    main(args)


if __name__ == '__main__':
    main_wrapper()
