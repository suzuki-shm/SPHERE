#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-27

import argparse
import numpy as np
import pandas as pd
from sphere.sphere_utils import load_depth_file
from sphere.sphere_utils import get_logger
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
    parser.add_argument("-m", "--model_type",
                        type=str,
                        nargs="?",
                        help="type of statistical model",
                        default="wrappedcauchy")
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
        pars = ["rho", "kappa"]
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
        density = linearcardioid_pdf(theta,
                                     loc=mu,
                                     rho=pars_values["rho"]["mean"])
    elif model == "cardioid":
        density = cardioid_pdf(theta,
                               loc=mu,
                               rho=pars_values["rho"]["mean"])
    elif model == "wrappedcauchy":
        density = wrappedcauchy(theta,
                                loc=mu,
                                rho=pars_values["rho"]["mean"])
    elif model == "vonmises":
        density = vonmises.pdf(theta,
                               loc=mu,
                               kappa=pars_values["kappa"]["mean"])
    elif model == "sslinearcardioid":
        density = sslinearcardioid_pdf(theta,
                                       loc=mu,
                                       rho=pars_values["rho"]["mean"],
                                       lambda_=pars_values["lambda"]["mean"])
    elif model == "sscardioid":
        density = sscardioid_pdf(theta,
                                 loc=mu,
                                 rho=pars_values["rho"]["mean"],
                                 lambda_=pars_values["lambda"]["mean"])
    elif model == "ssvonmises":
        density = ssvonmises_pdf(theta,
                                 loc=mu,
                                 kappa=pars_values["kappa"]["mean"],
                                 lambda_=pars_values["lambda"]["mean"])
    elif model == "sswrappedcauchy":
        density = sswrappedcauchy_pdf(theta,
                                      loc=mu,
                                      rho=pars_values["rho"]["mean"],
                                      lambda_=pars_values["lambda"]["mean"])
    result["mean"] = density

    # Min
    if model == "linearcardioid":
        density = linearcardioid_pdf(theta,
                                     loc=mu,
                                     rho=min(pars_values["rho"]["2.5%"],
                                             pars_values["rho"]["97.5%"]))
    elif model == "cardioid":
        density = cardioid_pdf(theta,
                               loc=mu,
                               rho=min(pars_values["rho"]["2.5%"],
                                       pars_values["rho"]["97.5%"]))
    elif model == "wrappedcauchy":
        density = wrappedcauchy(theta,
                                loc=mu,
                                rho=min(pars_values["rho"]["2.5%"],
                                        pars_values["rho"]["97.5%"]))
    elif model == "vonmises":
        density = vonmises.pdf(theta,
                               loc=mu,
                               kappa=min(pars_values["kappa"]["2.5%"],
                                         pars_values["kappa"]["97.5%"]))
    elif model == "sslinearcardioid":
        density = sslinearcardioid_pdf(theta,
                                       loc=mu,
                                       rho=min(pars_values["rho"]["2.5%"],
                                               pars_values["rho"]["97.5%"]),
                                       lambda_=min(pars_values["lambda"]["2.5%"],
                                                   pars_values["lambda"]["97.5%"]))
    elif model == "sscardioid":
        density = sscardioid_pdf(theta,
                                 loc=mu,
                                 rho=min(pars_values["rho"]["2.5%"],
                                         pars_values["rho"]["97.5%"]),
                                 lambda_=min(pars_values["lambda"]["2.5%"],
                                             pars_values["lambda"]["97.5%"]))
    elif model == "ssvonmises":
        density = ssvonmises_pdf(theta,
                                 loc=mu,
                                 kappa=min(pars_values["kappa"]["2.5%"],
                                           pars_values["kappa"]["97.5%"]),
                                 lambda_=min(pars_values["lambda"]["2.5%"],
                                             pars_values["lambda"]["97.5%"]))
    elif model == "sswrappedcauchy":
        density = sswrappedcauchy_pdf(theta,
                                      loc=mu,
                                      rho=min(pars_values["rho"]["2.5%"],
                                              pars_values["rho"]["97.5%"]),
                                      lambda_=min(pars_values["lambda"]["2.5%"],
                                                  pars_values["lambda"]["97.5%"]))
    result["min"] = density

    # Max
    if model == "linearcardioid":
        density = linearcardioid_pdf(theta,
                                     loc=mu,
                                     rho=max(pars_values["rho"]["2.5%"],
                                             pars_values["rho"]["97.5%"]))
    elif model == "cardioid":
        density = cardioid_pdf(theta,
                               loc=mu,
                               rho=max(pars_values["rho"]["2.5%"],
                                       pars_values["rho"]["97.5%"]))
    elif model == "wrappedcauchy":
        density = wrappedcauchy(theta,
                                loc=mu,
                                rho=max(pars_values["rho"]["2.5%"],
                                        pars_values["rho"]["97.5%"]))
    elif model == "vonmises":
        density = vonmises.pdf(theta,
                               loc=mu,
                               kappa=max(pars_values["kappa"]["2.5%"],
                                         pars_values["kappa"]["97.5%"]))
    elif model == "sslinearcarioid":
        density = sslinearcardioid_pdf(theta,
                                       loc=mu,
                                       rho=max(pars_values["rho"]["2.5%"],
                                               pars_values["rho"]["97.5%"]),
                                       lambda_=max(pars_values["lambda"]["2.5%"],
                                                   pars_values["lambda"]["97.5%"]))
    elif model == "sscardioid":
        density = sscardioid_pdf(theta,
                                 loc=mu,
                                 rho=max(pars_values["rho"]["2.5%"],
                                         pars_values["rho"]["97.5%"]),
                                 lambda_=max(pars_values["lambda"]["2.5%"],
                                             pars_values["lambda"]["97.5%"]))
    elif model == "ssvonmises":
        density = ssvonmises_pdf(theta,
                                 loc=mu,
                                 kappa=max(pars_values["kappa"]["2.5%"],
                                           pars_values["kappa"]["97.5%"]),
                                 lambda_=max(pars_values["lambda"]["2.5%"],
                                             pars_values["lambda"]["97.5%"]))
    elif model == "sswrappedcauchy":
        density = sswrappedcauchy_pdf(theta,
                                      loc=mu,
                                      rho=max(pars_values["rho"]["2.5%"],
                                              pars_values["rho"]["97.5%"]),
                                      lambda_=max(pars_values["lambda"]["2.5%"],
                                                  pars_values["lambda"]["97.5%"]))
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

    ax2._r_label_position._t = (22.5 + 180, 0.0)
    ax2._r_label_position.invalidate()

    for label in ax.get_yticklabels():
        ax.figure.texts.append(label)

    return ax2


def main(args, logger):
    fs = args["fs"]
    model = args["model_type"]
    index = args["index"]

    df = load_depth_file(args["depth_file_path"])
    summary_df = pd.read_csv(args["estimated_tsv"], sep="\t", index_col=0)
    length = len(df)
    x = np.linspace(-np.pi, np.pi, length)
    width = 2 * np.pi / length

    pars = get_target_parameter(model)
    mu_stats = get_mu_stats(summary_df)
    pars_stats = get_parameter_stats(summary_df, pars, index)
    density = get_density(model, pars_stats, mu_stats, length)

    fig = plt.figure(figsize=(10, 15))

    ax11 = fig.add_subplot(2, 1, 1)
    ax11.plot(x, density["mean"])
    ax11.fill_between(x, density["min"], density["max"], facecolor="pink")
    ax11.set_xlabel("Genomic position", fontsize=fs)
    ax11.set_ylabel("Probability density", fontsize=fs)
    ax11.tick_params(labelsize=fs)
    ax11.set_ylim(0, max(density["max"])*1.1)
    ax12 = ax11.twinx()
    ax12.bar(x, df["depth"], alpha=0.3, width=width)
    ax12.set_ylabel("Observed depth", fontsize=fs)
    ax12.tick_params(labelsize=fs)

    ax21 = fig.add_subplot(2, 1, 2, projection="polar")
    ax21.plot(x, density["mean"])
    ax21.fill_between(x, density["min"], density["max"], facecolor="pink")
    ax22 = polar_twin(ax21)
    ax22.bar(x, df["depth"], alpha=0.3, width=width)
    ax21.tick_params(labelsize=fs)
    ax21.set_theta_zero_location("N")
    ax22.tick_params(labelsize=fs)
    ax22.set_theta_zero_location("N")

    plt.savefig(args["output_dest"])


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
