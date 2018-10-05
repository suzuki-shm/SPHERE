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
from scipy.integrate import quad
from scipy.special import i0, i1, lpmv
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
                            "jonespewsey",
                            "dvonmises",
                            "sevonmises",
                            "sejonespewsey",
                            "invsevonmises",
                            "invsejonespewsey",
                            "miaecardioid",
                            "miaewrappedcauchy",
                            "miaevonmises",
                            "miaejonespewsey",
                            "invmiaecardioid",
                            "invmiaewrappedcauchy",
                            "invmiaevonmises",
                            "invmiaejonespewsey",
                            "invmievonmises"
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
    kappa_model = ("vonmises", "dvonmises", "jonespewsey", "djonespewsey")
    kappa_ae_model = ("miaevonmises", "miaejonespewsey",
                      "invmiaevonmises", "invmiaejonespewsey")
    kappa_se_model = ("sevonmises", "invsevonmises",
                      "sejonespewsey", "invsejonespewsey")
    jonespewsey = ("jonespewsey", "djonespewsey",
                   "sejonespewsey", "invsejonespewsey",
                   "miaejonespewsey", "invmiaejonespewsey")
    rho_model = ("linearcardioid", "cardioid", "wrappedcauchy")
    rho_ae_model = ("miaecardioid", "miaewrappedcauchy",
                    "invmiaecardioid", "invmiaewrappedcauchy")
    if model in kappa_model:
        pars = ["kappa"]
        if model in jonespewsey:
            pars.append("psi")
    elif model in kappa_ae_model:
        pars = ["kappa", "nu"]
        if model in jonespewsey:
            pars.append("psi")
    elif model in kappa_se_model:
        pars = ["kappa", "lambda"]
        if model in jonespewsey:
            pars.append("psi")
    elif model in rho_model:
        pars = ["rho"]
    elif model in rho_ae_model:
        pars = ["rho", "nu"]
    else:
        raise ValueError("Invalid input of model:{0}".format(model))

    return pars


def mix_density(density, alpha):
    d = (alpha * density).sum(axis=0)
    return d


def calc_tp(t, f, fd, theta):
    tp = t - f / fd
    return tp


# Basic distributions
def linearcardioid_pdf(theta, loc, rho):
    d = 1 / (2 * np.pi)
    d *= (1 + 2 * rho * (np.abs(np.abs(theta - loc) - np.pi) - np.pi / 2))
    return d


def cardioid_pdf(theta, loc, rho):
    d = 1 / (2 * np.pi) * (1 + 2 * rho * np.cos(theta - loc))
    return d


def wrappedcauchy_pdf(theta, loc, rho):
    d = (1 - np.power(rho, 2))
    m = 2 * np.pi * (1 +
                     np.power(rho, 2) -
                     2 * rho * np.cos(theta - loc))
    d = d / m
    return d


def jonespewsey_pdf(theta, loc, kappa, psi):
    def molecule(theta, loc, kappa, psi):
        d = np.power(np.cosh(kappa * psi) +
                     np.sinh(kappa * psi) * np.cos(theta - loc), 1/psi)
        return d
    if abs(psi) < 1e-10:
        p = vonmises.pdf(theta, kappa=kappa, loc=loc)
    else:
        m = molecule(theta, loc, kappa, psi)
        C = 2 * np.pi * lpmv(0, 1/psi, np.cosh(kappa*psi))
        p = m / C
    return p


# Papakonstanitinou transmation (or symmetric extended transformation)
def trans_se(theta, lambda_, loc):
    return theta-loc + lambda_ * np.sin(theta-loc)


def sevonmises_pdf(theta, loc, kappa, lambda_):
    def molecule(theta, loc, kappa, lambda_):
        return np.exp(kappa * np.cos(trans_se(theta, lambda_, loc)))
    p = molecule(theta, loc, kappa, lambda_)
    C = quad(molecule, -np.pi, np.pi, args=(0, kappa, lambda_))[0]
    return p / C


def sejonespewsey_pdf(theta, loc, kappa, psi, lambda_):
    def molecule(theta, loc, kappa, psi):
        d = np.power(np.cosh(kappa * psi) +
                     np.sinh(kappa * psi) *
                     np.cos(trans_se(theta, lambda_, loc)), 1/psi)
        return d
    p = molecule(theta, loc, kappa, lambda_)
    C = quad(molecule, -np.pi, np.pi, args=(0, kappa, lambda_))[0]
    return p / C


# inverse Batschelet transformation
# (or inverse symmetric extended transformation)
def trans_inv_se(theta, lambda_, loc):
    def inv_batschelet_trans_se(theta, lambda_, loc):
        t = theta
        for k in range(8):
            f = t - (1+lambda_) * np.sin(t-loc) / 2 - theta
            fd = 1 - (1+lambda_) * np.cos(t-loc) / 2
            tp = calc_tp(t, f, fd, theta)
            t = tp
        return t
    return ((1 - lambda_) / (1 + lambda_) * theta +
            2 * lambda_ / (1 + lambda_) *
            inv_batschelet_trans_se(theta, lambda_, loc))


def invsevonmises_pdf(theta, loc, kappa, lambda_):
    def normalized_constraint_inv_se(kappa, lambda_):
        def f(theta, kappa, lambda_):
            return ((1 - (1 + lambda_) * np.cos(theta) / 2) * vonmises.pdf((theta) - (1 - lambda_) * np.sin(theta) / 2, kappa=kappa, loc=0))
        return quad(f, -np.pi, np.pi, args=(kappa, lambda_))[0]

    # inverse transformation by Newton's method
    C = normalized_constraint_inv_se(kappa, lambda_)
    theta_trans = trans_inv_se(theta, lambda_, loc)
    p = vonmises.pdf(theta_trans, loc=loc, kappa=kappa) / C
    return p


def invsejonespewsey_pdf(theta, loc, kappa, psi, lambda_):
    def normalized_constraint_inv_se(kappa, psi, lambda_):
        def f(theta, kappa, psi, lambda_):
            if abs(psi) < 1e-10:
                return (1 - (1 + lambda_) * np.cos(theta) / 2) * vonmises.pdf((theta) - (1 - lambda_) * np.sin(theta) / 2, kappa=kappa, loc=0)
            else:
                return (1 - (1 + lambda_) * np.cos(theta) / 2) * jonespewsey_pdf((theta) - (1 - lambda_) * np.sin(theta) / 2, kappa=kappa, psi=psi, loc=0)
        return quad(f, -np.pi, np.pi, args=(kappa, psi, lambda_))[0]

    # inverse transformation by Newton's method
    C = normalized_constraint_inv_se(kappa, psi, lambda_)
    theta_trans = trans_inv_se(theta, lambda_, loc)
    p = jonespewsey_pdf(theta_trans, loc=loc, psi=psi, kappa=kappa) / C
    return p


# Mode invariance asymmetry extended transformation
def miaecardioid_pdf(theta, loc, rho, nu):
    theta_trans = theta - nu * np.sin(theta - loc) * np.sin(theta - loc)
    d = cardioid_pdf(theta_trans, loc=loc, rho=rho)
    return d


def miaewrappedcauchy_pdf(theta, loc, rho, nu):
    theta_trans = theta - nu * np.sin(theta - loc) * np.sin(theta - loc)
    d = wrappedcauchy_pdf(theta_trans, loc=loc, rho=rho)
    return d


def miaejonespewsey_pdf(theta, loc, kappa, nu, psi):
    theta_trans = theta - nu * np.sin(theta - loc) * np.sin(theta - loc)
    d = jonespewsey_pdf(theta_trans, loc=loc, kappa=kappa, psi=psi)
    return d


def miaevonmises_pdf(theta, loc, kappa, nu):
    theta_trans = theta - nu * np.sin(theta - loc) * np.sin(theta - loc)
    return vonmises.pdf(theta_trans, loc=loc, kappa=kappa)


# inverse mode invariance asymmetry extended transformation
def inv_trans_sin2(theta, loc, nu):
    # Inverse transformation by Newton's method
    t = theta
    for i in range(8):
        f = t + nu * np.sin(t - loc)**2 - theta
        fd = 2 * nu * np.sin(t - loc) * np.cos(t - loc)
        tp = calc_tp(t, f, fd, theta)
        t = tp
    return t


def invmiaecardioid_pdf(theta, loc, rho, nu):
    theta_trans = inv_trans_sin2(theta, loc, nu)
    d = cardioid_pdf(theta_trans, loc=loc, rho=rho)
    return d


def invmiaewrappedcauchy_pdf(theta, loc, rho, nu):
    theta_trans = inv_trans_sin2(theta, loc, nu)
    d = wrappedcauchy_pdf(theta_trans, loc=loc, rho=rho)
    return d


def invmiaejonespewsey_pdf(theta, loc, kappa, nu, psi):
    theta_trans = inv_trans_sin2(theta, loc, nu)
    d = jonespewsey_pdf(theta_trans, loc=loc, kappa=kappa, psi=psi)
    return d


def invmiaevonmises_pdf(theta, loc, kappa, nu):
    theta_trans = inv_trans_sin2(theta, loc, nu)
    return vonmises.pdf(theta_trans, loc=loc, kappa=kappa)


# Inverse Abe-Pewsey-Fujisawa transformation
def inv_trans_APF(theta, nu, lambda_, loc):
    t = theta
    for k in range(8):
        f = t-loc - nu * np.sin(t-loc) + lambda_ * np.power(np.sin(t-loc - nu * np.sin(t-loc)), 2)
        fd = (1 + 2 * lambda_ * np.sin(t-loc - nu * np.sin(t-loc)) * np.cos(t-loc - nu * np.sin(t-loc))) * (1 - nu * np.cos(t-loc))
        tp = calc_tp(t, f, fd, theta)
        t = tp
    return t


def invmievonmises_pdf(theta, kappa, nu, lambda_, loc):
    alpha1 = i1(kappa) / i0(kappa)
    # inverse transformation by Newton's method
    inv_theta = inv_trans_APF(theta, nu, lambda_, loc)
    C = (1 - nu * alpha1)
    p = vonmises.pdf(inv_theta, loc=0, kappa=kappa) / C
    return p


def get_density(model, pars_values, L, stat_type):
    theta = np.linspace(-np.pi, np.pi, L)

    mu = pars_values["mu"][stat_type]
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
    elif model == "jonespewsey":
        density = jonespewsey_pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type],
            psi=pars_values["psi"][stat_type]
        )
    elif model == "dvonmises":
        density = vonmises.pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type]
        )
        norm = np.linalg.norm(density, axis=1, keepdims=True, ord=1)
        density = density / norm
    elif model == "sevonmises":
        density = sevonmises_pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type],
            lambda_=pars_values["lambda"][stat_type]
        )
    elif model == "sejonespewsey":
        density = sejonespewsey_pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type],
            psi=pars_values["psi"][stat_type],
            lambda_=pars_values["lambda"][stat_type]
        )
    elif model == "invsevonmises":
        density = invsevonmises_pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type],
            lambda_=pars_values["lambda"][stat_type]
        )
    elif model == "invsejonespewsey":
        density = invsejonespewsey_pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type],
            psi=pars_values["psi"][stat_type],
            lambda_=pars_values["lambda"][stat_type]
        )
    elif model == "miaecardioid":
        density = miaecardioid_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"][stat_type],
            nu=pars_values["nu"][stat_type]
        )
    elif model == "miaewrappedcauchy":
        density = miaewrappedcauchy_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"][stat_type],
            nu=pars_values["nu"][stat_type]
        )
    elif model == "miaevonmises":
        density = miaevonmises_pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type],
            nu=pars_values["nu"][stat_type]
        )
    elif model == "miaejonespewsey":
        density = miaejonespewsey_pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type],
            psi=pars_values["psi"][stat_type],
            nu=pars_values["nu"][stat_type]
        )
    elif model == "invmiaecardioid":
        density = invmiaecardioid_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"][stat_type],
            nu=pars_values["nu"][stat_type]
        )
    elif model == "invmiaewrappedcauchy":
        density = invmiaewrappedcauchy_pdf(
            theta,
            loc=mu,
            rho=pars_values["rho"][stat_type],
            nu=pars_values["nu"][stat_type]
        )
    elif model == "invmiaevonmises":
        density = invmiaevonmises_pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type],
            nu=pars_values["nu"][stat_type]
        )
    elif model == "invmiaejonespewsey":
        density = invmiaejonespewsey_pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type],
            psi=pars_values["psi"][stat_type],
            nu=pars_values["nu"][stat_type]
        )
    elif model == "invmievonmises":
        density = invmiaevonmises_pdf(
            theta,
            loc=mu,
            kappa=pars_values["kappa"][stat_type],
            nu=pars_values["nu"][stat_type],
            lambda_=pars_values["lambda"][stat_type]
        )
    return density


def get_mu_stats(sdf, mode):
    if mode == "sampling":
        stats_type = ["2.5%", "mean", "97.5%"]
    elif mode == "optimizing":
        stats_type = ["mle"]

    pars_values = {}
    pars_values["mu"] = {}
    pars_df1 = sdf[sdf.index.str.match("O\[\d+,1\]")]
    pars_df2 = sdf[sdf.index.str.match("O\[\d+,2\]")]
    print(sdf)
    K = len(pars_df1)
    for st in stats_type:
        # The stats by Stan is not reliable, because ori parameter is not
        # linear value. So we have to calculate stats from O1 and O2
        O1 = pars_df1[st].values
        O2 = pars_df2[st].values
        v = np.arctan2(O1, O2)
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
    pars_df = sdf[sdf.index.str.match("alpha\[\d+\]")]
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
        # select parameter value raleted to input index
        pars_df = sdf[sdf.index.str.match("{0}\[{1},*\d*\]".format(p, index))]
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

    # Last index were discarded not to visualize overwrapped petal
    X_seg = np.linspace(0, 2 * np.pi, pn+1)[:-1]
    Y_seg = segment_depth(Y, pn)
    width = 2 * np.pi / (pn * 1.1)

    pars = get_target_parameter(model)
    mu_stats = get_mu_stats(sdf, mode)
    alpha_stats = get_alpha_stats(sdf, mode)
    pars_stats = get_parameter_stats(sdf, pars, i, mode)
    pars_stats.update(mu_stats)
    density = get_density(model,
                          pars_stats,
                          length,
                          mean_type)
    alpha = alpha_stats["alpha"][mean_type]
    weighted_density = alpha * density
    mean_density = mix_density(density, alpha)

    fig = plt.figure(figsize=(10, 15))

    ax11 = fig.add_subplot(2, 1, 1)
    ax11.tick_params(labelsize=fs)
    ax11.plot(X, mean_density)
    for d in weighted_density:
        ax11.plot(X, d)
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
    ax21.plot(X, mean_density)
    for d in weighted_density:
        ax21.plot(X, d)
    ax21.set_rticks(np.linspace(0, round(ax21.get_rmax()+0.05, 1), 3))
    ax21.set_theta_zero_location("N")
    ax22 = polar_twin(ax21)
    ax22.tick_params(labelsize=fs)
    ax22.bar(X_seg, Y_seg, alpha=0.3, width=width, align="edge")
    ax22.set_theta_zero_location("N")
    ax22.set_rticks(np.linspace(0, round(ax22.get_rmax(), 0), 3))
    ax22.set_rlabel_position(22.5 + 180)


def main(args, logger):
    fs = args["fs"]
    model = args["model_type"]
    index = args["index"]
    pn = args['pn']
    mode = args["M"]

    df = load_depth_file(args["depth_file_path"])
    summary_df = pd.read_csv(args["estimated_tsv"], sep="\t", index_col=0)
    plot_circular_dist(summary_df, df, fs, model, index, pn, mode)

    plt.savefig(args["output_dest"])


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
