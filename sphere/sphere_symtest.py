#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2018-04-26

from collections import namedtuple
from sphere.sphere_utils import load_depth_file
from sphere.sphere_cstats import (mean_resultant_length,
                                  mean_direction,
                                  sin_moment,
                                  cos_moment)
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


def perwey_test(theta, depth):
    """
    Computes the Perwey test on sample theta

    The Perwey test is a nonparametric test of the null hypothesis

    References
    ----------
    .. [1] Pewsey, A. "Testing Circular Symmetry". The Canadian Journal of
           Statistics. Vol. 30(2002): 591-600
    """
    PerweyResult = namedtuple("PerweyResult", ('statistic', 'pvalue'))
    n_depth = np.sum(depth)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    C = np.sum(depth * cos_theta)
    S = np.sum(depth * sin_theta)

    mrl = mean_resultant_length(C, S, np.sum(depth))
    md = mean_direction(S, C)
    b2 = sin_moment(depth, theta, p=2, loc=md)
    a2 = cos_moment(depth, theta, p=2, loc=md)
    a3 = cos_moment(depth, theta, p=3, loc=md)
    a4 = cos_moment(depth, theta, p=4, loc=md)
    var_b2 = ((1.0-a4)/2.0-2.0*a2+2.0*a2/mrl*(a3+(a2*(1.0-a2))/mrl))/n_depth

    z = b2 / np.sqrt(var_b2)
    p = 1 - norm.cdf(abs(z))
    p = 2 * np.min([p, 1-p])

    return PerweyResult(z, p)


def main(args):
    result = []
    for f in args["depth_file_path"]:
        file_name = os.path.basename(f)
        df = load_depth_file(f)
        length = len(df)
        x = np.arange(1, length+1, 1)
        theta = x / float(length) * 2.0 * np.pi
        depth = df["depth"].values
        z, p = perwey_test(theta, depth)

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
