#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-27

import argparse
import numpy as np
import pandas as pd
from sphere.sphere_utils import compress_depth, load_depth_file
from sphere.sphere_utils import get_logger
try:
    import matplotlib
    matplotlib.use("Agg")
finally:
    import matplotlib.pyplot as plt


def argument_parse(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("depth_file_path",
                        type=str,
                        help="file path of coverage depth")
    parser.add_argument("estimated_tsv",
                        type=str,
                        help="file path of estimated output tsv file")
    parser.add_argument("output_dest",
                        type=str,
                        help="destination of output figure")
    parser.add_argument("-cl", "--compressedlength",
                        dest="cl",
                        nargs="?",
                        default=10000,
                        type=int,
                        help="Compressed length of genome (default: 10000)")
    parser.add_argument("-fs", "--fontsize",
                        dest="fs",
                        type=int,
                        nargs="?",
                        default=18,
                        help="Font size of figure (default: 18)")
    args = parser.parse_args(argv)
    return vars(args)


def extend_array(y, w, L):
    o = np.ones((w, y.size))
    y_e = (o * y).T.ravel()[:L]
    return y_e


def main(args, logger):
    fs = args["fs"]
    df = load_depth_file(args["depth_file_path"])
    summary_df = pd.read_csv(args["estimated_tsv"], sep="\t", index_col=0)
    I = len(df)
    w = int(I / args["cl"])
    y = df["depth"]
    y_c = compress_depth(y, I, args["cl"])
    y_c_e = extend_array(y_c, w, I)
    x1 = np.arange(0, I, 1)
    x2 = np.arange(0, y_c.size, 1)

    m = "mean"
    l = "2.5%"
    u = "97.5%"
    y_eap = np.array([summary_df.loc["flex[{0}]".format(i), m] for i in x2])
    y_low = np.array([summary_df.loc["flex[{0}]".format(i), l] for i in x2])
    y_upp = np.array([summary_df.loc["flex[{0}]".format(i), u] for i in x2])
    t_eap = np.array([summary_df.loc["trend[{0}]".format(i), m] for i in x2])
    t_low = np.array([summary_df.loc["trend[{0}]".format(i), l] for i in x2])
    t_upp = np.array([summary_df.loc["trend[{0}]".format(i), u] for i in x2])
    l_eap = np.array([summary_df.loc["lambda[{0}]".format(i), m] for i in x2])
    l_low = np.array([summary_df.loc["lambda[{0}]".format(i), l] for i in x2])
    l_upp = np.array([summary_df.loc["lambda[{0}]".format(i), u] for i in x2])

    y_eap_e = extend_array(y_eap, w, I)
    y_low_e = extend_array(y_low, w, I)
    y_upp_e = extend_array(y_upp, w, I)
    t_eap_e = extend_array(t_eap, w, I)
    t_low_e = extend_array(t_low, w, I)
    t_upp_e = extend_array(t_upp, w, I)
    l_eap_e = extend_array(l_eap, w, I)
    l_low_e = extend_array(l_low, w, I)
    l_upp_e = extend_array(l_upp, w, I)

    fig = plt.figure(figsize=(10, 15))

    ax1 = fig.add_subplot(3, 1, 1)
    ax1.plot(x1, y_eap_e, label="estimated")
    ax1.fill_between(x1, y_low_e, y_upp_e, facecolor="pink")
    ax1.set_xlabel("Genomic position", fontsize=fs)
    ax1.set_ylabel("Cauchy trend", fontsize=fs)
    ax1.tick_params(labelsize=fs)
    ax1.legend(fontsize=fs)

    ax2 = fig.add_subplot(3, 1, 2)
    ax2.plot(x1, t_eap_e, label="estimated")
    ax2.fill_between(x1, t_low_e, t_upp_e, facecolor="pink")
    ax2.set_xlabel("Genomic position", fontsize=fs)
    ax2.set_ylabel("Triagonal trend", fontsize=fs)
    ax2.tick_params(labelsize=fs)
    ax2.legend(fontsize=fs)

    ax3 = fig.add_subplot(3, 1, 3)
    ax3.plot(x1, y, label="observed")
    ax3.plot(x1, y_c_e, label="filtered")
    ax3.plot(x1, l_eap_e, label="estimated")
    ax3.fill_between(x1, l_low_e, l_upp_e, facecolor="pink")
    ax3.set_xlabel("Genomic position", fontsize=fs)
    ax3.set_ylabel("Coverage depth", fontsize=fs)
    ax3.tick_params(labelsize=fs)
    ax3.legend(fontsize=fs)

    plt.savefig(args["output_dest"])


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
