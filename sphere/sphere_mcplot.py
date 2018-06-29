#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-27

import argparse
import numpy as np
import pandas as pd
from sphere.sphere_utils import load_depth_file
from sphere.sphere_utils import get_logger
from sphere.sphere_utils import segment_depth
from scipy.stats import vonmises
try:
    import matplotlib
    matplotlib.use("Agg")
finally:
    import matplotlib.pyplot as plt


def argument_parse(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dest",
                        type=str,
                        help="destination of output figure")
    parser.add_argument("depth_file_path",
                        type=str,
                        help="path of coverage depth")
    parser.add_argument("estimated_tsv",
                        type=str,
                        help="path of estimated output tsv file")
    parser.add_argument("index",
                        type=int,
                        help="index of visualized sample")
    parser.add_argument("-pn", "--petalnumber",
                        dest="pn",
                        type=int,
                        nargs="?",
                        default=50,
                        help="Petal number in rose diagram (default: 50)")
    parser.add_argument("-m", "--model_type",
                        nargs="?",
                        default="vonmises",
                        type=str,
                        choices=[
                            "linearcardioid",
                            "cardioid",
                            "wrappedcauchy",
                            "vonmises",
                            "sslinearcardioid",
                            "sscardioid",
                            "ssvonmises",
                            "sswrappedcauchy",
                            "statespacetrigonal",
                            "statespacelinear",
                            "trigonal",
                            "linear"
                        ],
                        help="type of statistical model",
                        )
    parser.add_argument("-M", "--method",
                        dest="M",
                        nargs="?",
                        default="sampling",
                        type=str,
                        choices=["sampling",
                                 "optimizing"],
                        help="adaptation method")
    parser.add_argument("-fs", "--fontsize",
                        dest="fs",
                        type=int,
                        nargs="?",
                        default=18,
                        help="Font size of figure (default: 18)")
    args = parser.parse_args(argv)
    return vars(args)


def get_target_parameter(model):
    if model == "linearcardioid":
        pars = ["rho"]
    elif model == "cardioid":
        pars = ["rho"]
    elif model == "wrappedcauchy":
        pars = ["rho"]
    elif model == "vonmises":
        pars = ["kappa"]
    elif model == "sslinearcardioid":
        pars = ["rho", "lambda"]
    elif model == "sscardioid":
        pars = ["rho", "lambda"]
    elif model == "sswrappedcauchy":
        pars = ["rho", "lambda"]
    elif model == "ssvonmises":
        pars = ["kappa", "lambda"]
    else:
        raise ValueError("Invalid input of model:{0}".format(model))

    return pars


def mix_density(density, alpha):
    d = (alpha * density).sum(axis=0)
    return d


def linearcardioid_pdf(theta, loc, rho):
    d = 1 / (2 * np.pi)
    d *= (1 + 2 * rho * (np.abs(np.abs(theta - loc) - np.pi) - np.pi / 2))
    return d


def cardioid_pdf(theta, loc, rho):
    d = 1 / (2 * np.pi) * (1 + 2 * rho * np.cos(theta - loc))
    return d


def wrappedcauchy_pdf(theta, loc, rho):
    d = (1 - np.power(rho, 2))
    m = (2 * np.pi * (1 + np.power(rho, 2) - 2 * rho * np.cos(theta - loc)))
    d = d / m
    return d


def sslinearcardioid_pdf(theta, loc, rho, lambda_):
    d = 1 / (2 * np.pi)
    d *= (1 + 2 * rho * (np.abs(np.abs(theta - loc) - np.pi) - np.pi / 2))
    d *= (1 + lambda_ * np.sin(theta - loc))
    return d


def sscardioid_pdf(theta, loc, rho, lambda_):
    d = 1 / (2 * np.pi) * (1 + 2 * rho * np.cos(theta - loc))
    d *= (1 + lambda_ * np.sin(theta - loc))
    return d


def sswrappedcauchy_pdf(theta, loc, rho, lambda_):
    d = (1 - np.power(rho, 2))
    m = (2 * np.pi * (1 + np.power(rho, 2) - 2 * rho * np.cos(theta - loc)))
    d = d / m
    d *= (1 + lambda_ * np.sin(theta - loc))
    return d


def ssvonmises_pdf(theta, loc, kappa, lambda_):
    d = vonmises.pdf(theta, loc=loc, kappa=kappa)
    d *= (1 + lambda_ * np.sin(theta - loc))
    return d


def get_density(model, pars_values, L, stat_type):
    theta = np.linspace(-np.pi, np.pi, L)

    mu = pars_values["mu"][stat_type]
    alpha = pars_values["alpha"][stat_type]
    # EAP
    if model == "linearcardioid":
        density = linearcardioid_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"][stat_type]
        )
    elif model == "cardioid":
        density = cardioid_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"][stat_type]
        )
    elif model == "wrappedcauchy":
        density = wrappedcauchy_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"][stat_type]
        )
    elif model == "vonmises":
        density = vonmises.pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type]
        )
    elif model == "sslinearcardioid":
        density = sslinearcardioid_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"][stat_type],
            lambda_=pars_values["lambda"][stat_type]
        )
    elif model == "sscardioid":
        density = sscardioid_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"][stat_type],
            lambda_=pars_values["lambda"][stat_type]
        )
    elif model == "ssvonmises":
        density = ssvonmises_pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type],
            lambda_=pars_values["lambda"][stat_type]
        )
    elif model == "sswrappedcauchy":
        density = sswrappedcauchy_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"][stat_type],
            lambda_=pars_values["lambda"][stat_type]
        )
    density = mix_density(density, alpha)

    return density


def get_mu_stats(sdf, mode):
    if mode == "sampling":
        stats_type = ["2.5%", "mean", "97.5%"]
    elif mode == "optimizing":
        stats_type = ["mle"]

    pars_values = {}
    pars_values["mu"] = {}
    pars_df1 = sdf[sdf.index.str.match("O\[.+,0\]")]
    pars_df2 = sdf[sdf.index.str.match("O\[.+,1\]")]
    K = len(pars_df1)
    for st in stats_type:
        # The stats by Stan is not reliable, because ori parameter is not
        # linear value. So we have to calculate stats from O1 and O2
        O1 = pars_df1[st].values
        O2 = pars_df2[st].values
        # The range of this value is -pi to pi.
        # When visualizing the coverage, we'll use 0 to 2pi.
        # To transform the value, add pi
        v = np.arctan2(O1, O2) + np.pi
        v = v.reshape(K, 1)
        pars_values["mu"][st] = v
    return pars_values


def get_alpha_stats(sdf, mode):
    if mode == "sampling":
        stats_type = ["2.5%", "mean", "97.5%"]
    elif mode == "optimizing":
        stats_type = ["mle"]

    pars_values = {}
    pars_values["alpha"] = {}
    pars_df = sdf[sdf.index.str.match("alpha\[.+\]")]
    K = len(pars_df)
    for st in stats_type:
        v = pars_df[st].values
        v = v.reshape(K, 1)
        pars_values["alpha"][st] = v
    return pars_values


def get_parameter_stats(sdf, pars, index, mode):
    if mode == "sampling":
        stats_type = ["2.5%", "mean", "97.5%"]
    elif mode == "optimizing":
        stats_type = ["mle"]

    pars_values = {}
    for p in pars:
        pars_values[p] = {}
        pars_df = sdf[sdf.index.str.match("{0}\[{1},.+\]".format(p, index))]
        K = len(pars_df)
        for st in stats_type:
            v = pars_df[st].values
            v = v.reshape(K, 1)
            pars_values[p][st] = v
    return pars_values


def polar_twin(ax):
    ax2 = ax.figure.add_axes(ax.get_position(),
                             projection='polar',
                             label='twin',
                             frameon=False,
                             theta_direction=ax.get_theta_direction(),
                             theta_offset=ax.get_theta_offset())
    ax2.xaxis.set_visible(False)

    for label in ax.get_yticklabels():
        ax.figure.texts.append(label)

    return ax2


def plot_circular_dist(sdf, depth_df, fs, model, i, pn, mode):
    length = len(depth_df)
    X = np.linspace(0, 2*np.pi, length)
    Y = depth_df["depth"]
    if mode == "sampling":
        mean_type = "mean"
        min_type = "2.5%"
        max_type = "97.5%"
    elif mode == "optimizing":
        mean_type = "mle"
    xaxis_range = np.linspace(1, length, 5)
    xaxis_label = ["{:.1e}".format(l) for l in xaxis_range]

    X_seg = np.linspace(0, 2 * np.pi, pn)
    Y_seg = segment_depth(Y, pn)
    width = 2 * np.pi / (pn * 1.1)

    pars = get_target_parameter(model)
    mu_stats = get_mu_stats(sdf, mode)
    alpha_stats = get_alpha_stats(sdf, mode)
    pars_stats = get_parameter_stats(sdf, pars, i, mode)
    pars_stats.update(mu_stats)
    pars_stats.update(alpha_stats)
    mean_density = get_density(model,
                               pars_stats,
                               length,
                               mean_type)

    fig = plt.figure(figsize=(10, 15))

    ax11 = fig.add_subplot(2, 1, 1)
    ax11.tick_params(labelsize=fs)
    ax11.plot(X, mean_density, color="#ed7d31")
    ax11.set_xlabel("Genomic position", fontsize=fs)
    ax11.set_ylabel("Probability density", fontsize=fs)
    ax11.set_ylim(bottom=0)
    ax11.set_xticks(np.linspace(0, 2*np.pi, 5))
    ax11.set_xticklabels(xaxis_label)
    ax12 = ax11.twinx()
    ax12.tick_params(labelsize=fs)
    ax12.plot(X, Y, alpha=0.3)
    ax12.set_ylim(bottom=0)
    ax12.set_ylabel("Observed depth", fontsize=fs)

    ax21 = fig.add_subplot(2, 1, 2, projection="polar")
    ax21.tick_params(labelsize=fs)
    ax21.plot(X, mean_density, color="#ed7d31")
    ax21.set_rticks(np.linspace(0, round(ax21.get_rmax()+0.05, 1), 3))
    ax21.set_theta_zero_location("N")
    ax22 = polar_twin(ax21)
    ax22.tick_params(labelsize=fs)
    ax22.bar(X_seg, Y_seg, alpha=0.3, width=width)
    ax22.set_theta_zero_location("N")
    ax22.set_rticks(np.linspace(0, round(ax22.get_rmax(), 0), 3))
    ax22.set_rlabel_position(22.5 + 180)


def plot_statespace(sdf, depth_df, fs, model, i, pn, mode):
    length = len(depth_df)
    X = depth_df["location"] - 1
    Y = depth_df["depth"]
    if mode == "sampling":
        m = "mean"
        low = "2.5%"
        high = "97.5%"
    elif mode == "optimizing":
        m = "mle"
    T = np.linspace(0, 2*np.pi, length)
    width = 2 * np.pi / length
    xaxis_range = np.linspace(1, length, 5)
    xaxis_label = ["{:.1e}".format(l) for l in xaxis_range]

    t_m = np.array([sdf.loc["trend[{0},{1}]".format(i, x), m] for x in X])
    l_m = np.array([sdf.loc["lambda[{0},{1}]".format(i, x), m] for x in X])
    if mode == "sampling":
        t_l = [sdf.loc["trend[{0},{1}]".format(i, x), low] for x in X]
        l_l = [sdf.loc["lambda[{0},{1}]".format(i, x), low] for x in X]
        t_h = [sdf.loc["trend[{0},{1}]".format(i, x), high] for x in X]
        l_h = [sdf.loc["lambda[{0},{1}]".format(i, x), high] for x in X]
        t_l = np.array(t_l)
        l_l = np.array(l_l)
        t_h = np.array(t_h)
        l_h = np.array(l_h)

    fig = plt.figure(figsize=(10, 20))

    ax1 = fig.add_subplot(3, 1, 1)
    ax1.tick_params(labelsize=fs)
    ax1.plot(X, t_m, label=m, color="#ed7d31")
    ax1.fill_between(X, t_l, t_h, facecolor="#ff9e00", alpha="0.3")
    ax1.set_xlabel("Genomic position", fontsize=fs)
    ax1.set_ylabel("Trend", fontsize=fs)
    ax1.set_xticks(xaxis_range)
    ax1.set_xticklabels(xaxis_label)

    ax21 = fig.add_subplot(3, 1, 2)
    ax21.plot(X, l_m, label=m, color="#ed7d31")
    if mode == "sampling":
        ax21.fill_between(X, l_l, l_h, facecolor="#ff9e00", alpha=0.3)
    ax21.tick_params(labelsize=fs)
    ax21.set_ylabel("Potential", fontsize=fs)
    ax21.set_ylim(bottom=0)
    ax22 = ax21.twinx()
    ax22.tick_params(labelsize=fs)
    ax22.plot(X, Y, label="observed")
    ax22.set_xlabel("Genomic position", fontsize=fs)
    ax22.set_ylabel("Coverage depth", fontsize=fs)
    ax22.set_xticks(xaxis_range)
    ax22.set_xticklabels(xaxis_label)
    ax22.set_ylim(bottom=0)

    ax31 = fig.add_subplot(3, 1, 3, projection="polar")
    ax31.tick_params(labelsize=fs)
    ax31.plot(T, l_m, color="#ed7d31")
    if mode == "sampling":
        ax31.fill_between(T, l_l, l_h, facecolor="#ff9e00", alpha=0.3)
    ax31.set_rticks(np.linspace(0, round(ax31.get_rmax()+0.05, 1), 3))
    ax31.set_theta_zero_location("N")
    ax32 = polar_twin(ax31)
    ax32.tick_params(labelsize=fs)
    ax32.bar(T, Y, alpha=0.3, width=width)
    ax32.set_theta_zero_location("N")
    ax32.set_rticks(np.linspace(0, round(ax32.get_rmax(), 0), 3))
    ax32.set_rlabel_position(22.5 + 180)


def main(args, logger):
    fs = args["fs"]
    model = args["model_type"]
    index = args["index"]
    pn = args['pn']
    mode = args["M"]

    df = load_depth_file(args["depth_file_path"])
    summary_df = pd.read_csv(args["estimated_tsv"], sep="\t", index_col=0)
    circular_model = ["linearcardioid",
                      "cardioid",
                      "wrappedcauchy",
                      "vonmises",
                      "sslinearcardioid",
                      "sscardioid",
                      "ssvonmises",
                      "sswrappedcauchy"]
    statespace_model = ["statespacetrigonal",
                        "statespacelinear",
                        "trigonal",
                        "linear"]
    if model in circular_model:
        plot_circular_dist(summary_df, df, fs, model, index, pn, mode)
    elif model in statespace_model:
        plot_statespace(summary_df, df, fs, model, index, pn, mode)

    plt.savefig(args["output_dest"])


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
