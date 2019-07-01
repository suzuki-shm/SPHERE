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
    if len(df) == 0:
        msg = "File {0} is empty".format(depth_file_path)
        raise ValueError(msg)
    return df


def load_multiple_depth_file(depth_file_path: list):
    list_ = []
    for i, f in enumerate(depth_file_path):
        df = load_depth_file(f)
        df["subject"] = i+1
        list_.append(df)
    c_df = pd.concat(list_)
    return c_df


def moving_filter(d: pd.Series, s: int=None, w: int=None, ftype="median") -> pd.Series:
    if w != 1:
        # Append head and foot part of array, considering circular structure
        d_head = d[:np.floor(w/2).astype(int)]
        d_foot = d[-np.floor(w/2).astype(int):]
        dr = d_foot.append(d).append(d_head)
        # Take rolling median
        if ftype == "median":
            dr = dr.rolling(window=w, min_periods=1).median()
        elif ftype == "sum":
            dr = dr.rolling(window=w, min_periods=1).sum()
        elif ftype == "mvariance":
            dr_median = dr.rolling(window=w, min_periods=1).median()
            dr_std = dr.rolling(window=w, min_periods=1).std()
            dr[(dr_median - dr).abs() > dr_std] = np.nan
        else:
            raise ValueError("Invalid filter type: {0}".format(ftype))
        # Drop out double calculated parts
        chi = w-1
        if w % 2 != 0:
            dr = dr[chi:]
        else:
            dr = dr[chi:-1]
    else:
        dr = d
    # Take results with stride length
    dr = dr.reset_index(drop=True)
    dr = dr[list(range(0, dr.size, s))]
    dr = dr.reset_index(drop=True)
    dr = dr.round()
    try:
        dr = dr.astype(int)
    except ValueError:
        logger = get_logger()
        logger.warning("Compressed depth still have NaN")
    return dr


def compress_length(dl: int, s: int, w: int) -> int:
    cl = np.ceil(dl / s)
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
                "PTR",  "wPTR", "mwPTR", "MRL", "CV", "CSD"]
    elif model_name == "explinearcardioid":
        pars = ["alpha", "O", "rho", "ori",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "cardioid":
        pars = ["alpha", "O", "rho", "kappa", "ori",
                "PTR",  "wPTR", "mwPTR", "MRL", "CV", "CSD"]
    elif model_name == "vonmises" or model_name == "vonmisesloss":
        pars = ["alpha", "O", "kappa", "ori",
                "PTR",  "wPTR", "mwPTR", "MRL", "CV", "CSD"]
    elif model_name == "wrappedcauchy":
        pars = ["alpha", "O", "rho", "kappa", "ori",
                "PTR",  "wPTR", "mwPTR", "MRL", "CV", "CSD"]
    elif model_name == "jonespewsey":
        pars = ["alpha", "O", "kappa", "psi", "ori",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "dvonmises":
        pars = ["alpha", "O", "kappa", "ori",
                "PTR",  "wPTR", "mwPTR", "MRL", "CV", "CSD"]
    elif model_name == "sevonmises":
        pars = ["alpha", "O", "kappa", "ori", "lambda",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "sejonespewsey":
        pars = ["alpha", "O", "kappa", "psi", "ori", "lambda",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "invsevonmises":
        pars = ["alpha", "O", "kappa", "ori", "lambda",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "invsejonespewsey":
        pars = ["alpha", "O", "kappa", "psi", "ori", "lambda",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "miaecardioid":
        pars = ["alpha", "O", "rho", "kappa", "ori", "nu",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "miaevonmises":
        pars = ["alpha", "O", "kappa", "ori", "nu",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "miaewrappedcauchy":
        pars = ["alpha", "O", "rho", "kappa", "ori", "nu",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "miaejonespewsey":
        pars = ["alpha", "O", "kappa", "psi", "ori", "nu",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "invmiaecardioid":
        pars = ["alpha", "O", "rho", "kappa", "ori", "nu",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "invmiaevonmises":
        pars = ["alpha", "O", "kappa", "ori", "nu",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "invmiaewrappedcauchy":
        pars = ["alpha", "O", "rho", "kappa", "ori", "nu",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "invmiaejonespewsey":
        pars = ["alpha", "O", "kappa", "psi", "ori", "nu",
                "PTR",  "wPTR", "mwPTR"]
    elif model_name == "invmievonmises":
        pars = ["alpha", "O", "kappa", "ori", "nu", "lambda",
                "PTR",  "wPTR", "mwPTR"]

    if has_log_lik:
        pars.append("log_lik")

    return pars


def main():
    pass


if __name__ == '__main__':
    main()
