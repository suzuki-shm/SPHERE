#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

import pandas as pd
import numpy as np


def load_depth_file(depth_file_path):
    df = pd.read_csv(depth_file_path,
                     sep="\t",
                     index_col=0,
                     names=["position", "depth"])
    return df


def compress_depth(v, I, cl, t="ceil"):
    w = int(I / cl)
    if w == 0:
        w = 1

    if t == "ceil":
        cl_comp = np.ceil(I/w).astype(int)
    elif t == "floor":
        cl_comp = np.floor(I/w).astype(int)
    v_resized = np.resize(v, (cl_comp, w))
    v_median = np.median(v_resized, axis=1)
    v_compressed = np.round(v_median).astype(int)
    return v_compressed


def main():
    pass


if __name__ == '__main__':
    main()
