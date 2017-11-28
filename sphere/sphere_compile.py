#!/usr/bin/env python
# vim:fileencoding=utf-8
# Author: Shinya Suzuki
# Created: 2017-09-26

import argparse
import pystan
from sphere.stan_utils import save_model
from sphere.sphere_utils import get_logger


def argument_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_path",
                        type=str,
                        help="file path of compiled stan model")
    parser.add_argument("--model_type", "-m",
                        type=str,
                        choices=["trigonal",
                                 "linear",
                                 "mix_vonmises",
                                 "vonmises"],
                        default="trigonal",
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
                int S ;
                int L ;
                int<lower=1, upper=L> LOCATION[I] ;
                int<lower=1, upper=S> SUBJECT[I] ;
                int<lower=0> DEPTH[I] ;
            }

            parameters {
                unit_vector[2] O ;
                real<lower=0> H[S] ;
                real flex0[S] ;
                vector<lower=-pi()/2, upper=pi()/2>[L-1] flex_raw[S] ;
                real<lower=0> sigma_flex[S] ;
            }

            transformed parameters{
                real<lower=-2*pi(), upper=2*pi()> ori ;
                vector[L] flex[S] ;
                vector[L] trend[S] ;
                vector<lower=0>[L] lambda[S] ;

                // convert unit vector
                ori = atan2(O[1], O[2]) ;

                for(s in 1:S){
                    // flex
                    flex[s, 1] = flex0[s] ;
                    for(l in 2:L){
                        flex[s, l] = flex[s, l-1] + sigma_flex[s] * tan(flex_raw[s, l-1]) ;
                    }

                    // trend from replication rate
                    for(l in 1:L){
                        trend[s, l] = H[s] / 2.0 * (cos(l * 2.0 * pi() / L - ori) + 1.0) ;
                    }
                    lambda[s] = exp(flex[s] + trend[s]) ;
                }
            }

            model {
                for(i in 1:I){
                    DEPTH[i] ~ poisson(lambda[SUBJECT[i], LOCATION[i]]) ;
                }
            }

            generated quantities {
                real<lower=1.0> PTR[S] ;
                vector[I] log_lik ;

                for(s in 1:S){
                    PTR[s] = exp(H[s]) ;
                }
                for(i in 1:I){
                    log_lik[i] = poisson_lpmf(DEPTH[i] | lambda[SUBJECT[i], LOCATION[i]]) ;
                }
            }
        """
    elif model == "linear":
        model_code = """
            data {
                int I ;
                int S ;
                int L ;
                int<lower=1, upper=L> LOCATION[I] ;
                int<lower=1, upper=S> SUBJECT[I] ;
                int<lower=0> DEPTH[I] ;
            }

            parameters {
                unit_vector[2] O ;
                real<lower=0> H[S] ;
                real flex0[S] ;
                vector<lower=-pi()/2, upper=pi()/2>[I-1] flex_raw[S] ;
                real<lower=0> sigma_flex[S] ;
            }

            transformed parameters{
                real<lower=-2*pi(), upper=2*pi()> ori ;
                vector[L] flex[S] ;
                vector[L] trend[S] ;
                vector<lower=0>[L] lambda[S] ;

                // convert unit vector
                ori = atan2(O[1], O[2]) ;

                for(s in 1:S){
                    // flex
                    flex[s, 1] = flex0[s] ;
                    for(l in 2:L){
                        flex[s, l] = flex[s, l-1] + sigma_flex[s] * tan(flex_raw[s, l-1]) ;
                    }

                    // trend from replication rate
                    for(l in 1:L){
                        trend[s, l] = 2.0 * H[s] / I * fabs(fabs(l * 2.0 * pi() / L - ori ) - L / 2.0) ;
                    }
                    lambda[s] = exp(flex[s] + trend[s]) ;
                }
            }

            model {
                for(i in 1:I){
                    DEPTH[i] ~ poisson(lambda[SUBJECT[i], LOCATION[i]]) ;
                }
            }

            generated quantities {
                real<lower=1.0> PTR[S] ;
                vector[I] log_lik ;

                for(s in 1:S){
                    PTR[s] = exp(H[s]) ;
                }
                for(i in 1:I){
                    log_lik[i] = poisson_lpmf(DEPTH[i] | lambda[SUBJECT[i], LOCATION[i]]) ;
                }
            }
        """
    elif model == "vonmises":
        model_code = """
            data {
                int I ;
                int S ;
                int L ;
                int<lower=1, upper=L> LOCATION[I] ;
                int<lower=1, upper=S> SUBJECT[I] ;
                int<lower=0> DEPTH[I] ;
            }

            transformed data {
                real RADIAN[I] ;
                for (i in 1:I){
                    RADIAN[i] = 2.0 * pi() * LOCATION[i] / L ;
                }
            }

            parameters {
                unit_vector[2] O ;
                real<lower=0> kappa[S] ;
            }

            transformed parameters{
                real<lower=-2*pi(), upper=2*pi()> ori ;

                // convert unit vector
                ori = atan2(O[1], O[2]) ;
            }

            model {
                for(i in 1:I){
                    target += DEPTH[i] * von_mises_lpdf(RADIAN[i] | ori, kappa[SUBJECT[i]]) ;
                }
            }

            generated quantities {
                real MRL[S] ;
                real CV[S] ;
                real CSD[S] ;
                vector[I] log_lik ;

                for(s in 1:S){
                    MRL[s] = modified_bessel_first_kind(1, kappa[s]) / modified_bessel_first_kind(0, kappa[s]) ;
                    CV[s] = 1 - MRL[s] ;
                    CSD[s] = sqrt(-2 * log(MRL[s])) ;
                }
                for(i in 1:I){
                    log_lik[i] = DEPTH[i] * von_mises_lpdf(RADIAN[i] | ori, kappa[SUBJECT[i]]) ;
                }
            }
        """
    model = pystan.StanModel(model_code=model_code)
    if output_path is not None:
        save_model(output_path, model)
    return model


def main(args, logger):
    args = argument_parse()
    compile_model(args["output_path"], args["model_type"])
    logger.info("Stan model is compiled to {0}.".format(args["output_path"]))


def main_wrapper():
    args = argument_parse()
    logger = get_logger(__name__)
    main(args, logger)


if __name__ == '__main__':
    main_wrapper()
