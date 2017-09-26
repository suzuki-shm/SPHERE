#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

import argparse
import pystan
import pickle
from logging import getLogger, DEBUG, Formatter, StreamHandler


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
    parser.add_argument("model_path",
                        type=str,
                        help="file path of stan model")
    parser.add_argument("output_path",
                        type=str,
                        help="file path of compiled stan model")
    parser.set_defaults(trans=False)
    args = parser.parse_args()
    return vars(args)


def compile_model(model_path, output_path):
    # Stanのモデルを読み込んでコンパイルする
    model = pystan.StanModel(file=model_path)
    with open(output_path, "wb") as f:
        pickle.dump(model, f)
    return output_path


def main(args, logger):
    output_path = compile_model(args["model_path"], args["output_path"])
    logger.info("Stan model is compiled to {0}.".format(output_path))


if __name__ == '__main__':
    main(argument_parse(), get_logger())
