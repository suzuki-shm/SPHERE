#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

import pandas as pd
import numpy as np
from logging import getLogger, INFO, Formatter, StreamHandler


def load_depth_file(depth_file_path: str):
    df = pd.read_csv(depth_file_path,
                     sep="\t",
                     names=["genome", "location", "depth"])
    if df["genome"].unique().size != 1:
        raise ValueError("File contains multiple mapping result")
    x = np.arange(1, len(df)+1, 1)
    f_df = pd.DataFrame(x, columns=["location"])
    j_df = df.merge(f_df, on=["location"], how="outer")
    genome_name = df["genome"].unique()[0]
    j_df["depth"] = j_df["depth"].fillna(0)
    j_df["depth"] = j_df["depth"].astype(int)
    j_df["genome"] = genome_name
    return j_df


def load_multiple_depth_file(depth_file_path: list):
    list_ = []
    for i, f in enumerate(depth_file_path):
        df = load_depth_file(f)
        df["subject"] = i+1
        df = df[["subject", "location", "depth"]]
        list_.append(df)
    c_df = pd.concat(list_)
    return c_df


def compress_depth(d: pd.Series, s: int=None, w: int=None) -> pd.Series:
    if w != 1:
        # Append head and foot part of array, considering circular structure
        d_head = d[:(w-1)]
        d_foot = d[-(w-1):]
        dr = d_foot.append(d).append(d_head)
        # Take rolling median
        dr = dr.rolling(window=w, min_periods=w).median()
        dr = dr.dropna().reset_index(drop=True)
        # Drop out double calculated parts
        chi = np.ceil((w-1)/2).astype(int)
        cfi = np.floor((w-1)/2).astype(int)
        dr = dr[chi:-cfi].reset_index(drop=True)
    else:
        dr = d
    dr = dr[list(range(0, dr.size, s))].reset_index(drop=True)
    dr = dr.round().astype(int)
    return dr


def compress_length(dl: int, s: int, w: int) -> int:
    if w != 1:
        fl = dl - w % 2
    else:
        fl = dl
    cl = np.ceil(fl / s)
    return cl


def segment_depth(v: np.ndarray, cl: int) -> np.ndarray:
    length = v.size
    w1 = window_length(length, cl)
    w2 = w1 + 1

    A = np.array([[w1, w2], [1, 1]])
    B = np.array([v.size, cl])
    m, n = np.round(np.linalg.solve(A, B)).astype(int)

    v1 = v[:w1 * m]
    v1_resized = np.resize(v1, (m, w1))
    v1_median = np.median(v1_resized, axis=1)
    v1_compressed = np.round(v1_median).astype(int)
    if n != 0:
        v2 = v[w2 * n:]
        v2_resized = np.resize(v2, (n, w2))
        v2_median = np.median(v2_resized, axis=1)
        v2_compressed = np.round(v2_median).astype(int)
        vr = np.r_[v1_compressed, v2_compressed]
    else:
        vr = v1_compressed
    return vr


def window_length(I, cl):
    w = int(I / cl)
    if w == 0:
        w = 1
    return w


def get_logger(name=None):
    if name is None:
        logger = getLogger(__name__)
    else:
        logger = getLogger(name)
    logger.setLevel(INFO)
    log_fmt = '%(asctime)s : %(name)s : %(levelname)s : %(message)s'
    formatter = Formatter(log_fmt)
    stream_handler = StreamHandler()
    stream_handler.setLevel(INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    return logger


def get_pars(model_name, has_log_lik=False):
    if model_name == "linearcardioid":
        pars = ["alpha", "O", "rho", "ori",
                "PTR",  "mPTR", "wmPTR", "MRL", "CV", "CSD"]
    elif model_name == "cardioid":
        pars = ["alpha", "O", "rho", "kappa", "ori",
                "PTR",  "mPTR", "wmPTR", "MRL", "CV", "CSD"]
    elif model_name == "vonmises":
        pars = ["alpha", "O", "kappa", "ori",
                "PTR",  "mPTR", "wmPTR", "MRL", "CV", "CSD"]
    elif model_name == "wrappedcauchy":
        pars = ["alpha", "O", "rho", "kappa", "ori",
                "PTR",  "mPTR", "wmPTR", "MRL", "CV", "CSD"]
    elif model_name == "jonespewsey":
        pars = ["alpha", "O", "kappa", "psi", "ori",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "dvonmises":
        pars = ["alpha", "O", "kappa", "ori",
                "PTR",  "mPTR", "wmPTR", "MRL", "CV", "CSD"]
    elif model_name == "sevonmises":
        pars = ["alpha", "O", "kappa", "ori", "lambda",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "sejonespewsey":
        pars = ["alpha", "O", "kappa", "psi", "ori", "lambda",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "invsevonmises":
        pars = ["alpha", "O", "kappa", "ori", "lambda",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "invsejonespewsey":
        pars = ["alpha", "O", "kappa", "psi", "ori", "lambda",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "miaecardioid":
        pars = ["alpha", "O", "rho", "kappa", "ori", "nu",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "miaevonmises":
        pars = ["alpha", "O", "kappa", "ori", "nu",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "miaewrappedcauchy":
        pars = ["alpha", "O", "rho", "kappa", "ori", "nu",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "miaejonespewsey":
        pars = ["alpha", "O", "kappa", "psi", "ori", "nu",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "invmiaecardioid":
        pars = ["alpha", "O", "rho", "kappa", "ori", "nu",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "invmiaevonmises":
        pars = ["alpha", "O", "kappa", "ori", "nu",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "invmiaewrappedcauchy":
        pars = ["alpha", "O", "rho", "kappa", "ori", "nu",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "invmiaejonespewsey":
        pars = ["alpha", "O", "kappa", "psi", "ori", "nu",
                "PTR",  "mPTR", "wmPTR"]
    elif model_name == "invmievonmises":
        pars = ["alpha", "O", "kappa", "ori", "nu", "lambda",
                "PTR",  "mPTR", "wmPTR"]

    if has_log_lik:
        pars.append("log_lik")

    return pars


def main():
    pass


if __name__ == '__main__':
    main()
