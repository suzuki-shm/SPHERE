#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

import argparse
import numpy as np
from sphere.sphere_utils import load_depth_file
from sphere.sphere_utils import segment_depth
from sphere.sphere_utils import get_logger
from matplotlib import gridspec
try:
    import matplotlib
    matplotlib.use('Agg')
finally:
    import matplotlib.pyplot as plt


def argument_parse(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("depth_file_path",
                        type=str,
                        help="file path of coverage depth")
    parser.add_argument("output_dest",
                        type=str,
                        help="destination of output file")
    parser.add_argument("-pn", "--petalnumber",
                        dest="pn",
                        type=int,
                        nargs="?",
                        default=50,
                        help="Petal number in rose diagram (default: 50)")
    parser.add_argument("-fs", "--fontsize",
                        dest="fs",
                        type=int,
                        nargs="?",
                        default=30,
                        help="Font size of figure (default: 18)")
    args = parser.parse_args(argv)
    return vars(args)


def main(args, logger):
    fs = args["fs"]
    df = load_depth_file(args["depth_file_path"])
    length = len(df)
    X = np.arange(0, length, 1)
    Y = df["depth"].values
    pn = args["pn"]

    t1 = np.arange(0, 2*np.pi, 2*np.pi/pn)
    t2 = np.arange(0, 2*np.pi, 2*np.pi/length)
    Y_seg = segment_depth(Y, pn)
    width = 2*np.pi / (pn+10)

    fig = plt.figure(figsize=(20, 20))
    gs = gridspec.GridSpec(2, 2)

    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(X, Y)
    ax1.set_xlabel("Genomic position", fontsize=fs)
    ax1.set_ylabel("Coverage depth", fontsize=fs)
    ax1.tick_params(labelsize=fs)

    ax2 = fig.add_subplot(gs[1, 0], projection="polar")
    ax2.bar(t1, Y_seg, width=width, align="edge")
    ax2.set_theta_zero_location("N")
    ax2.set_xticks(np.arange(0, 360, 360/6) / 360 * 2 * np.pi)
    ax2.set_xticklabels(np.arange(0, length, length/6, dtype=int))
    ax2.tick_params(labelsize=fs)

    ax3 = fig.add_subplot(gs[1, 1], projection="polar")
    ax3.plot(t2, Y)
    ax3.set_theta_zero_location("N")
    ax3.set_xticks(np.arange(0, 360, 360/6) / 360 * 2 * np.pi)
    ax3.set_xticklabels(np.arange(0, length, length/6, dtype=int))
    ax3.tick_params(labelsize=fs)

    plt.savefig(args["output_dest"])


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
