#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

from logging import getLogger, DEBUG, Formatter, StreamHandler
import os
import argparse
import pandas as pd
import numpy as np
import pickle


def get_logger():
    logger = getLogger(__name__)
    logger.setLevel(DEBUG)
    log_fmt = '%(asctime)s : %(name)s : %(levelname)s : %(message)s'
    formatter = Formatter(log_fmt)
    stream_handler = StreamHandler()
    stream_handler.setLevel(DEBUG)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    return logger


def argument_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("depth_file_path",
                        type=str,
                        help="file path of coverage depth")
    parser.add_argument("compiled_model_path",
                        type=str,
                        help="file path of compiled Stan model")
    parser.add_argument("output_path",
                        type=str,
                        help="file path of output tsv file")
    parser.add_argument("-fop", "--fit_out_path",
                        nargs="?",
                        default=None,
                        type=str,
                        help="file path of estimated result (default: None)")
    parser.add_argument("-cl", "--compressedlength",
                        nargs="?",
                        default=10000,
                        type=int,
                        help="Compressed length of genome (default: 10000)")
    parser.add_argument("-si", "--staniter",
                        nargs="?",
                        default=3000,
                        type=int,
                        help="Number of Stan iteration (default: 3000)")
    parser.add_argument("-sw", "--stanwarmup",
                        nargs="?",
                        default=1000,
                        type=int,
                        help="Number of Stan warm up (defaultL 1000)")
    parser.add_argument("-sc", "--stanchain",
                        nargs="?",
                        default=3,
                        type=int,
                        help="Number of Stan chain (default: 3)")
    parser.add_argument("-st", "--stanthin",
                        nargs="?",
                        default=1,
                        type=int,
                        help="Number of Stan thin (default: 1)")
    parser.add_argument("-ss", "--stanseed",
                        nargs="?",
                        default=1234,
                        type=int,
                        help="Number of Stan seed (defaultt: 1234)")
    parser.add_argument("-ff",
                        dest="ff",
                        action="store_true",
                        help="overwrite output file (default: False)")
    parser.set_defaults(trans=False)
    args = parser.parse_args()
    return vars(args)


def load_depth_file(depth_file_path):
    df = pd.read_csv(depth_file_path,
                     sep="\t",
                     index_col=0,
                     names=["position", "depth"])
    return df


def load_model(compiled_model_path):
    with open(compiled_model_path, "rb") as f:
        model = pickle.load(f)
    return model


def compress_depth(v, I, cl):
    w = int(I / cl)
    if w == 0:
        w = 1
    cl_ceil = np.ceil(I/w).astype(int)
    v_resized = np.resize(v, (cl_ceil, w))
    v_median = np.median(v_resized, axis=1)
    v_compressed = np.round(v_median).astype(int)
    return v_compressed


def sampling(model, v_c, I, si, sw, sc, st, ss):
    stan_data = {
        "I": I,
        "D": v_c
    }
    fit = model.sampling(data=stan_data,
                         iter=si,
                         warmup=sw,
                         chain=sc,
                         thin=st,
                         seed=ss)
    return fit


def summarize_fit(fit):
    summary = fit.summary()
    summary_df = pd.DataFrame(summary["summary"],
                              index=summary["summary_rownames"],
                              columns=summary["summary_colnames"])
    return summary_df


def check_output_path(ff, op):
    if ff is False:
        if os.path.exist(op):
            raise ValueError("{0} has been already existed. ".format(op),
                             "Use -ff to overwrite.")
        else:
            return 0
    else:
        if os.path.exist(op):
            return 1
        else:
            return 0


def save_fit(fit, fop):
    with open(fop, "wb") as f:
        pickle.dump(fit, f)


def main(args, logger):
    ope = check_output_path(args["ff"], args["output_path"])
    if ope:
        logger.warning("{0} will be overwrited".format(args["output_path"]))
    df = load_depth_file(args["depth_file_path"])
    model = load_model(args["compiled_model_path"])
    I = len(df)
    v_c = compress_depth(df["depth"], I, args["cl"])
    fit = sampling(model,
                   v_c,
                   I,
                   args["si"], args["sw"], args["sc"], args["st"], args["ss"])
    sdf = summarize_fit(fit)
    sdf.to_csv(args["output_path"], sep="\t")
    if args["fop"] is not None:
        save_fit(fit, args["fop"])


if __name__ == '__main__':
    main(argument_parse(), get_logger())
