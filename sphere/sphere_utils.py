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


def compress_depth(v: np.ndarray, I: int, cl: int, t="ceil"):
    w = window_length(I, cl)
    if t == "ceil":
        cl_comp = np.ceil(I / w).astype(int)
    elif t == "floor":
        cl_comp = np.floor(I / w).astype(int)
    v_resized = np.resize(v, (cl_comp, w))
    v_median = np.median(v_resized, axis=1)
    v_compressed = np.round(v_median).astype(int)
    return v_compressed


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
