#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

import argparse
import pystan
from logging import getLogger, DEBUG, Formatter, StreamHandler
from sphere.stan_utils import save_model


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
    parser.add_argument("output_path",
                        type=str,
                        help="file path of compiled stan model")
    parser.set_defaults(trans=False)
    args = parser.parse_args()
    return vars(args)


def compile_model(output_path=None, model="trigonal"):
    # Stanのモデルを読み込んでコンパイルする
    if model == "trigonal":
        model_code = """
            data {
                int I ;
                int D[I] ;
            }

            parameters {
                real flex0 ;
                real<lower=0> H ;
                real<lower=-1, upper=1> O[2] ;
                vector<lower=-pi()/2, upper=pi()/2>[I-1] flex_raw ;
                real<lower=0> sigma_flex ;
            }

            transformed parameters{
                vector<lower=0>[I] lambda ;
                vector[I] flex ;
                vector[I] trend ;

                // flex
                flex[1] = flex0 ;
                for(i in 2:I){
                    flex[i] = flex[i-1] + sigma_flex * tan(flex_raw[i-1]) ;
                }

                // trend from replication rate
                for(i in 1:I){
                    trend[i] = H / 2.0 * (cos((2.0 * pi() * i) / I - atan2(O[1], O[2])) + 1.0) ;
                }
                lambda = exp(flex + trend) ;

            }

            model {
                D ~ poisson(lambda) ;
            }

            generated quantities {
                real<lower=1.0> PTR ;
                vector[I] log_lik ;

                PTR = exp(H) ;
                for(i in 1:I){
                    log_lik[i] = poisson_lpmf(D[i] | lambda[i]) ;
                }
            }
        """
    elif model == "linear":
        model_code = """
            data {
                int I ;
                int D[I] ;
            }

            parameters {
                real flex0 ;
                real<lower=0> H ;
                real<lower=-1, upper=1> O[2] ;
                vector<lower=-pi()/2, upper=pi()/2>[I-1] flex_raw ;
                real<lower=0> sigma_flex ;
            }

            transformed parameters{
                vector<lower=0>[I] lambda ;
                vector[I] flex ;
                vector[I] trend ;

                // flex
                flex[1] = flex0 ;
                for(i in 2:I){
                    flex[i] = flex[i-1] + sigma_flex * tan(flex_raw[i-1]) ;
                }

                // trend from replication rate
                for(i in 1:I){
                    trend[i] = 2.0 * H / I * fabs(fabs(i - atan2(O[1], O[2]) / 2.0 / pi() * I) - I / 2.0) ;
                }
                lambda = exp(flex + trend) ;

            }

            model {
                D ~ poisson(lambda) ;
            }

            generated quantities {
                real<lower=1.0> PTR ;
                vector[I] log_lik ;

                PTR = exp(H) ;
                for(i in 1:I){
                    log_lik[i] = poisson_lpmf(D[i] | lambda[i]) ;
                }
            }
        """
    elif model == "mix":
        model_code = """
            data {
                int I ;
                int D[I] ;
            }

            parameters {
                real flex0 ;
                real<lower=0> H ;
                real<lower=-1, upper=1> O[2] ;
                vector<lower=-pi()/2, upper=pi()/2>[I-1] flex_raw ;
                real<lower=0> sigma_flex ;
                real<lower=0> p ;
            }

            transformed parameters{
                vector<lower=0>[I] lambda ;
                vector[I] flex ;
                vector[I] trigonal;
                vector[I] linear;
                vector[I] trend ;

                // flex
                flex[1] = flex0 ;
                for(i in 2:I){
                    flex[i] = flex[i-1] + sigma_flex * tan(flex_raw[i-1]) ;
                }

                // trend from replication rate
                for(i in 1:I){
                    trigonal[i] = (cos((2.0 * pi() * i) / I - atan2(O[1], O[2])) + 1.0) / 2.0 ;
                    linear[i] = 2.0 / I * fabs(fabs(i - atan2(O[1], O[2]) / 2.0 / pi() * I) - I / 2.0) ;
                    trend[i] = H * trigonal[i] * linear[i] ;
                }

                lambda = exp(flex + trend) ;
            }

            model {
                D ~ poisson(lambda) ;
            }

            generated quantities {
                real<lower=1.0> PTR ;
                vector[I] log_lik ;

                PTR = exp(H) ;
                for(i in 1:I){
                    log_lik[i] = poisson_lpmf(D[i] | lambda[i]) ;
                }
            }
        """
    model = pystan.StanModel(model_code=model_code)
    if output_path is not None:
        save_model(output_path, model)
    return model


def main(args, logger):
    args = argument_parse()
    compile_model(args["output_path"])
    logger.info("Stan model is compiled to {0}.".format(args["output_path"]))


def main_wrapper():
    args = argument_parse()
    logger = get_logger()
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
