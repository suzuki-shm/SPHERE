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
    parser.add_argument("depth_file_path",
                        type=str,
                        help="path of coverage depth")
    parser.add_argument("estimated_tsv",
                        type=str,
                        help="path of estimated output tsv file")
    parser.add_argument("output_dest",
                        type=str,
                        help="destination of output figure")
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


def linearcardioid_pdf(theta, loc, rho):
    d = 1 / (2 * np.pi)
    d *= (1 + 2 * rho * (np.abs(np.abs(theta - loc) - np.pi) - np.pi / 2))
    return d


def cardioid_pdf(theta, loc, rho):
    d = 1 / (2 * np.pi) * (1 + 2 * rho * np.cos(theta - loc))
    return d


def wrappedcauchy(theta, loc, rho):
    d = (1 - np.power(rho, 2))
    d /= (2 * np.pi * (1 + np.power(rho, 2) - 2 * rho * np.cos(theta - loc)))
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
    d /= (2 * np.pi * (1 + np.power(rho, 2) - 2 * rho * np.cos(theta - loc)))
    d *= (1 + lambda_ * np.sin(theta - loc))
    return d


def ssvonmises_pdf(theta, loc, kappa, lambda_):
    d = vonmises.pdf(theta, loc=loc, kappa=kappa)
    d *= (1 + lambda_ * np.sin(theta - loc))
    return d


def get_density(model, pars_values, mu_values, I):
    theta = np.linspace(-np.pi, np.pi, I)
    result = {}
    mu = mu_values["mu"]["mean"]

    # EAP
    if model == "linearcardioid":
        density = linearcardioid_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"]["mean"]
        )
    elif model == "cardioid":
        density = cardioid_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"]["mean"]
        )
    elif model == "wrappedcauchy":
        density = wrappedcauchy(
            theta,
            loc=mu,
            rho=pars_values["rho"]["mean"]
        )
    elif model == "vonmises":
        density = vonmises.pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"]["mean"]
        )
    elif model == "sslinearcardioid":
        density = sslinearcardioid_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"]["mean"],
            lambda_=pars_values["lambda"]["mean"]
        )
    elif model == "sscardioid":
        density = sscardioid_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"]["mean"],
            lambda_=pars_values["lambda"]["mean"]
        )
    elif model == "ssvonmises":
        density = ssvonmises_pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"]["mean"],
            lambda_=pars_values["lambda"]["mean"]
        )
    elif model == "sswrappedcauchy":
        density = sswrappedcauchy_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"]["mean"],
            lambda_=pars_values["lambda"]["mean"]
        )
    result["mean"] = density

    # Min
    if model == "linearcardioid":
        density = linearcardioid_pdf(
            theta,
            loc=mu,
            rho=min(pars_values["rho"]["2.5%"],
                    pars_values["rho"]["97.5%"])
        )
    elif model == "cardioid":
        density = cardioid_pdf(
            theta,
            loc=mu,
            rho=min(pars_values["rho"]["2.5%"],
                    pars_values["rho"]["97.5%"])
        )
    elif model == "wrappedcauchy":
        density = wrappedcauchy(
            theta,
            loc=mu,
            rho=min(pars_values["rho"]["2.5%"],
                    pars_values["rho"]["97.5%"])
        )
    elif model == "vonmises":
        density = vonmises.pdf(
            theta,
            loc=mu,
            kappa=min(pars_values["kappa"]["2.5%"],
                      pars_values["kappa"]["97.5%"])
        )
    elif model == "sslinearcardioid":
        density = sslinearcardioid_pdf(
            theta,
            loc=mu,
            rho=min(pars_values["rho"]["2.5%"],
                    pars_values["rho"]["97.5%"]),
            lambda_=min(pars_values["lambda"]["2.5%"],
                        pars_values["lambda"]["97.5%"])
        )
    elif model == "sscardioid":
        density = sscardioid_pdf(
            theta,
            loc=mu,
            rho=min(pars_values["rho"]["2.5%"],
                    pars_values["rho"]["97.5%"]),
            lambda_=min(pars_values["lambda"]["2.5%"],
                        pars_values["lambda"]["97.5%"])
        )
    elif model == "ssvonmises":
        density = ssvonmises_pdf(
            theta,
            loc=mu,
            kappa=min(pars_values["kappa"]["2.5%"],
                      pars_values["kappa"]["97.5%"]),
            lambda_=min(pars_values["lambda"]["2.5%"],
                        pars_values["lambda"]["97.5%"])
        )
    elif model == "sswrappedcauchy":
        density = sswrappedcauchy_pdf(
            theta,
            loc=mu,
            rho=min(pars_values["rho"]["2.5%"],
                    pars_values["rho"]["97.5%"]),
            lambda_=min(pars_values["lambda"]["2.5%"],
                        pars_values["lambda"]["97.5%"])
        )
    result["min"] = density

    # Max
    if model == "linearcardioid":
        density = linearcardioid_pdf(
            theta,
            loc=mu,
            rho=max(pars_values["rho"]["2.5%"],
                    pars_values["rho"]["97.5%"])
        )
    elif model == "cardioid":
        density = cardioid_pdf(
            theta,
            loc=mu,
            rho=max(pars_values["rho"]["2.5%"],
                    pars_values["rho"]["97.5%"])
        )
    elif model == "wrappedcauchy":
        density = wrappedcauchy(
            theta,
            loc=mu,
            rho=max(pars_values["rho"]["2.5%"],
                    pars_values["rho"]["97.5%"])
        )
    elif model == "vonmises":
        density = vonmises.pdf(
            theta,
            loc=mu,
            kappa=max(pars_values["kappa"]["2.5%"],
                      pars_values["kappa"]["97.5%"])
        )
    elif model == "sslinearcarioid":
        density = sslinearcardioid_pdf(
            theta,
            loc=mu,
            rho=max(pars_values["rho"]["2.5%"],
                    pars_values["rho"]["97.5%"]),
            lambda_=max(pars_values["lambda"]["2.5%"],
                        pars_values["lambda"]["97.5%"])
        )
    elif model == "sscardioid":
        density = sscardioid_pdf(
            theta,
            loc=mu,
            rho=max(pars_values["rho"]["2.5%"],
                    pars_values["rho"]["97.5%"]),
            lambda_=max(pars_values["lambda"]["2.5%"],
                        pars_values["lambda"]["97.5%"])
        )
    elif model == "ssvonmises":
        density = ssvonmises_pdf(
            theta,
            loc=mu,
            kappa=max(pars_values["kappa"]["2.5%"],
                      pars_values["kappa"]["97.5%"]),
            lambda_=max(pars_values["lambda"]["2.5%"],
                        pars_values["lambda"]["97.5%"])
        )
    elif model == "sswrappedcauchy":
        density = sswrappedcauchy_pdf(
            theta,
            loc=mu,
            rho=max(pars_values["rho"]["2.5%"],
                    pars_values["rho"]["97.5%"]),
            lambda_=max(pars_values["lambda"]["2.5%"],
                        pars_values["lambda"]["97.5%"])
        )
    result["max"] = density
    return result


def get_mu_stats(summary_df):
    stats_type = ["2.5%", "mean", "97.5%"]
    pars_values = {}
    pars_values["mu"] = {}
    for st in stats_type:
        O1 = summary_df.loc["O[0]", st]
        O2 = summary_df.loc["O[1]", st]
        mu = np.arctan2(O1, O2) + np.pi
        pars_values["mu"][st] = mu
    return pars_values


def get_parameter_stats(summary_df, pars, index):
    stats_type = ["2.5%", "mean", "97.5%"]
    pars_values = {}
    for p in pars:
        pars_values[p] = {}
    for st in stats_type:
        for p in pars:
            v = summary_df.loc["{0}[{1}]".format(p, index), st]
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


def plot_circular_dist(sdf, depth_df, fs, model, i, pn):
    length = len(depth_df)
    X = np.linspace(0, 2*np.pi, length)
    Y = depth_df["depth"]
    xaxis_range = np.linspace(1, length, 5)
    xaxis_label = ["{:.1e}".format(l) for l in xaxis_range]

    X_seg = np.linspace(-np.pi, np.pi, pn)
    Y_seg = segment_depth(Y, pn)
    width = 2 * np.pi / (pn*1.1)

    pars = get_target_parameter(model)
    mu_stats = get_mu_stats(sdf)
    pars_stats = get_parameter_stats(sdf, pars, i)
    density = get_density(model, pars_stats, mu_stats, length)

    fig = plt.figure(figsize=(10, 15))

    ax11 = fig.add_subplot(2, 1, 1)
    ax11.tick_params(labelsize=fs)
    ax11.plot(X, density["mean"], color="#ed7d31")
    ax11.fill_between(X,
                      density["min"],
                      density["max"],
                      facecolor="#ff9e00", alpha=0.3)
    ax11.set_xlabel("Genomic position", fontsize=fs)
    ax11.set_ylabel("Probability density", fontsize=fs)
    ax11.set_ylim(bottom=0, top=max(density["max"])*1.1)
    ax11.set_xticks(np.linspace(-np.pi, np.pi, 5))
    ax11.set_xticklabels(xaxis_label)
    ax12 = ax11.twinx()
    ax12.tick_params(labelsize=fs)
    ax12.bar(X_seg, Y_seg, alpha=0.3, width=width)
    ax12.set_ylim(bottom=0)
    ax12.set_ylabel("Observed depth", fontsize=fs)

    ax21 = fig.add_subplot(2, 1, 2, projection="polar")
    ax21.tick_params(labelsize=fs)
    ax21.plot(X, density["mean"], color="#ed7d31")
    ax21.fill_between(X,
                      density["min"],
                      density["max"],
                      facecolor="#ff9e00", alpha=0.3)
    ax21.set_rticks(np.linspace(0, round(ax21.get_rmax()+0.05, 1), 3))
    ax21.set_theta_zero_location("N")
    ax22 = polar_twin(ax21)
    ax22.tick_params(labelsize=fs)
    ax22.bar(X_seg, Y_seg, alpha=0.3, width=width)
    ax22.set_theta_zero_location("N")
    ax22.set_rticks(np.linspace(0, round(ax22.get_rmax(), 0), 3))
    ax22.set_rlabel_position(22.5 + 180)


def plot_statespace(sdf, depth_df, fs, model, i, pn):
    length = len(depth_df)
    X = np.arange(0, length, 1)
    Y = depth_df["depth"]
    T = np.linspace(0, 2*np.pi, length)
    width = 2 * np.pi / length
    xaxis_range = np.linspace(1, length, 5)
    xaxis_label = ["{:.1e}".format(l) for l in xaxis_range]

    m = "mean"
    low = "2.5%"
    high = "97.5%"
    t_eap = np.array([sdf.loc["trend[{0},{1}]".format(i, x), m] for x in X])
    t_l = np.array([sdf.loc["trend[{0},{1}]".format(i, x), low] for x in X])
    t_h = np.array([sdf.loc["trend[{0},{1}]".format(i, x), high] for x in X])
    l_eap = np.array([sdf.loc["lambda[{0},{1}]".format(i, x), m] for x in X])
    l_l = np.array([sdf.loc["lambda[{0},{1}]".format(i, x), low] for x in X])
    l_h = np.array([sdf.loc["lambda[{0},{1}]".format(i, x), high] for x in X])

    fig = plt.figure(figsize=(10, 20))

    ax1 = fig.add_subplot(3, 1, 1)
    ax1.tick_params(labelsize=fs)
    ax1.plot(X, t_eap, label="EAP", color="#ed7d31")
    ax1.fill_between(X, t_l, t_h, facecolor="#ff9e00", alpha="0.3")
    ax1.set_xlabel("Genomic position", fontsize=fs)
    ax1.set_ylabel("Trend", fontsize=fs)
    ax1.set_xticks(xaxis_range)
    ax1.set_xticklabels(xaxis_label)

    ax21 = fig.add_subplot(3, 1, 2)
    ax21.plot(X, l_eap, label="EAP", color="#ed7d31")
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
    ax31.plot(T, l_eap, color="#ed7d31")
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
        plot_circular_dist(summary_df, df, fs, model, index, pn)
    elif model in statespace_model:
        plot_statespace(summary_df, df, fs, model, index, pn)

    plt.savefig(args["output_dest"])


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
