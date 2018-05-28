#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-28

from sphere.sphere_utils import load_depth_file
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


def mean_resultant_length(C, S, n):
    # n means number of sample
    rho = np.sqrt(C**2 + S**2) / n
    return rho


def circular_variance(rho):
    V = 1 - rho
    return V


def circular_standard_deviation(rho):
    v = np.sqrt(-2.0 * np.log(rho))
    return v


def mean_direction(S, C):
    if C > 0 and S >= 0:
        md = np.arctan(S / C)
    elif C == 0 and S > 0:
        md = np.pi / 2.0
    elif C < 0:
        md = np.arctan(S/C) + np.pi
    elif C == 0 and S < 0:
        md = 3.0 * np.pi / 2.0
    else:
        md = np.arctan(S / C) + 2 * np.pi
    return md


def mean_direction_1dimension(md, I):
    md_1d = md / (2.0 * np.pi) * float(I)
    if md_1d < 0:
        return md_1d + float(I)
    else:
        return md_1d


def sin_moment(d, t, p=1, loc=0):
    Sbar = np.sum(d * np.sin(p * (t - loc))) / np.sum(d)
    return Sbar


def cos_moment(d, t, p=1, loc=0):
    Cbar = np.sum(d * np.cos(p * (t - loc))) / np.sum(d)
    return Cbar


def circular_skewness(d, t, mrl):
    beta2bar = sin_moment(d, t, p=2, loc=mrl)
    cs = beta2bar / np.power(1-mrl, 1.5)
    return cs


def circular_kurtosis(d, t, mrl):
    alpha2bar = cos_moment(d, t, p=2, loc=mrl)
    ck = (alpha2bar-np.power(mrl, 4)) / np.power(1-mrl, 2)
    return ck


def main(args):
    result = []
    for f in args["depth_file_path"]:
        file_name = os.path.basename(f)
        df = load_depth_file(f)
        length = len(df)
        x = np.arange(1, length+1, 1)
        theta = x / float(length) * 2.0 * np.pi
        depth = df["depth"].values
        n_depth = np.sum(depth)
        C = np.sum(depth * np.cos(theta))
        S = np.sum(depth * np.sin(theta))

        mrl = mean_resultant_length(C, S, n_depth)
        cv = circular_variance(mrl)
        csd = circular_standard_deviation(mrl)
        md = mean_direction(S, C)
        md_1d = mean_direction_1dimension(md, length)
        cs = circular_skewness(depth, theta, mrl)
        ck = circular_kurtosis(depth, theta, mrl)

        tmp = {
            "filepath": f,
            "filename": file_name,
            "mean_resultant_length": mrl,
            "circular_variance": cv,
            "circular_standard_deviation": csd,
            "mean_direction": md,
            "mean_position": md_1d,
            "circular_skewness": cs,
            "circular_kurtosis": ck
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
                            "circular_skewness",
                            "circular_kurtosis",
                            "filepath"
                         ]]
    result_df.to_csv(args["output_dest"], sep="\t")


def main_wrapper():
    args = argument_parse()
    main(args)


if __name__ == '__main__':
    main_wrapper()
