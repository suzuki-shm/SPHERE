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
                     names=["genome", "location", "depth"])
    if df["genome"].unique().size != 1:
        raise ValueError("File contains multiple mapping result")
    x = np.arange(1, len(df)+1, 1)
    f_df = pd.DataFrame(x, columns=["location"])
    j_df = df.merge(f_df, on=["location"], how="outer")
    genome_name = df["genome"].unique()[0]
    j_df["depth"] = j_df["depth"].fillna(0)
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


def compress_depth(v: np.ndarray, cl: int):
    nrows = v.size - cl + 1
    v1 = np.arange(nrows)[:, None]
    v2 = np.arange(cl)
    v_strided = v[v1 + v2]
    v_median = np.median(v_strided, axis=0)
    return np.round(v_median).astype(int)


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
