#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

import argparse
import numpy as np
from logging import getLogger, DEBUG, Formatter, StreamHandler
from sphere.sphere_utils import load_depth_file
from sphere.sphere_utils import compress_depth
from matplotlib import gridspec
try:
    import matplotlib
    matplotlib.use('Agg')
finally:
    import matplotlib.pyplot as plt


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


def argument_parse(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("depth_file_path",
                        type=str,
                        help="file path of coverage depth")
    parser.add_argument("output_dest",
                        type=str,
                        help="destination of output file")
    parser.add_argument("-cl", "--compressedlength",
                        dest="cl",
                        nargs="?",
                        default=10000,
                        type=int,
                        help="Compressed length of genome (default: 10000)")
    parser.add_argument("-np", "--npetal",
                        dest="np",
                        type=int,
                        nargs="?",
                        default=30,
                        help="Number of petal in rose diagram (default: 30)")
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
    I = len(df)
    x = df.index.values
    y = df["depth"].values
    y_f = compress_depth(y, I, args["np"], "floor")
    y_c = compress_depth(y, I, args["cl"])
    x_c = np.arange(1, y_c.size+1, 1)
    t1 = np.arange(0, 2*np.pi, 2*np.pi/args["np"])
    t2 = np.arange(0, 2*np.pi, 2*np.pi/I)
    width = 2 * np.pi / (args["np"]+10)

    fig = plt.figure(figsize=(20, 30))
    gs = gridspec.GridSpec(3, 2)

    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(x, y)
    ax1.set_xlabel("Genomic position", fontsize=fs)
    ax1.set_ylabel("Coverage depth", fontsize=fs)
    ax1.tick_params(labelsize=fs)

    ax2 = fig.add_subplot(gs[1, :])
    ax2.plot(x_c, y_c)
    ax2.set_xlabel("Compressed genomic position", fontsize=fs)
    ax2.set_ylabel("Compressed coverage depth", fontsize=fs)
    ax2.tick_params(labelsize=fs)

    ax3 = fig.add_subplot(gs[2, 0], projection="polar")
    ax3.bar(t1, y_f, width=width)
    ax3.set_theta_zero_location("N")
    ax3.set_xticks(np.arange(0, 360, 360/6) / 360 * 2 * np.pi)
    ax3.set_xticklabels(np.arange(0, I, I/6, dtype=int))
    ax3.tick_params(labelsize=fs)

    ax4 = fig.add_subplot(gs[2, 1], projection="polar")
    ax4.plot(t2, y)
    ax4.set_theta_zero_location("N")
    ax4.set_xticks(np.arange(0, 360, 360/6) / 360 * 2 * np.pi)
    ax4.set_xticklabels(np.arange(0, I, I/6, dtype=int))
    ax4.tick_params(labelsize=fs)

    plt.savefig(args["output_dest"])


def main_wrapper():
    args = argument_parse()
    logger = get_logger()
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
