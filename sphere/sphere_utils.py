#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

import pandas as pd
import numpy as np
from logging import getLogger, DEBUG, Formatter, StreamHandler


def load_depth_file(depth_file_path: str):
    df = pd.read_csv(depth_file_path,
                     sep="\t",
                     names=["genome", "position", "depth"])
    if df["genome"].unique().size != 1:
        raise ValueError("File contains multiple mapping result")
    x = np.arange(1, len(df)+1, 1)
    f_df = pd.DataFrame(x, columns=["position"])
    j_df = df.merge(f_df, on=["position"], how="outer")
    genome_name = df["genome"].unique()[0]
    j_df["depth"] = j_df["depth"].fillna(0)
    j_df["genome"] = genome_name
    return j_df


def compress_depth(v: np.ndarray, cl: int):
    I = v.size
    w1 = window_length(I, cl)
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


def get_logger(name):
    logger = getLogger(name)
    logger.setLevel(DEBUG)
    log_fmt = '%(asctime)s : %(name)s : %(levelname)s : %(message)s'
    formatter = Formatter(log_fmt)
    stream_handler = StreamHandler()
    stream_handler.setLevel(DEBUG)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    return logger


def main():
    pass


if __name__ == '__main__':
    main()
