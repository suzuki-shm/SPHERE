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
                int D[I] ;
            }

            parameters {
                real flex0 ;
                real<lower=0> H ;
                unit_vector[2] O ;
                vector<lower=-pi()/2, upper=pi()/2>[I-1] flex_raw ;
                real<lower=0> sigma_flex ;
            }

            transformed parameters{
                vector<lower=0>[I] lambda ;
                vector[I] flex ;
                vector[I] trend ;
                real<lower=-pi(), upper=pi()> ori ;

                // convert unit vector
                ori = atan2(O[1], O[2]) ;

                // flex
                flex[1] = flex0 ;
                for(i in 2:I){
                    flex[i] = flex[i-1] + sigma_flex * tan(flex_raw[i-1]) ;
                }

                // trend from replication rate
                for(i in 1:I){
                    trend[i] = H / 2.0 * (cos(i * 2.0 * pi() / I - ori) + 1.0) ;
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
                unit_vector[2] O ;
                vector<lower=-pi()/2, upper=pi()/2>[I-1] flex_raw ;
                real<lower=0> sigma_flex ;
            }

            transformed parameters{
                vector<lower=0>[I] lambda ;
                vector[I] flex ;
                vector[I] trend ;
                real<lower=-pi(), upper=pi()> ori ;

                // convert unit vector
                ori = atan2(O[1], O[2]) ;

                // flex
                flex[1] = flex0 ;
                for(i in 2:I){
                    flex[i] = flex[i-1] + sigma_flex * tan(flex_raw[i-1]) ;
                }

                // trend from replication rate
                for(i in 1:I){
                    trend[i] = 2.0 * H / I * fabs(fabs(i * 2.0 * pi() / I - ori ) - I / 2.0) ;
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
    elif model == "mix_vonmises":
        model_code = """
            data {
                int I ;
                int D[I] ;
            }

            transformed data {
                real R[I] ;
                for (i in 1:I){
                    R[i] = 2.0 * pi() * i / I ;
                }
            }

            parameters {
                simplex[3] theta ;
                unit_vector[2] O ;
                real<lower=0> kappa1 ;
                real<lower=0> kappa2 ;
                real<lower=0> kappa3 ;
            }

            transformed parameters{
                real<lower=-pi(), upper=pi()> ori ;

                // convert unit vector
                ori = atan2(O[1], O[2]) ;
            }

            model {
                real ps[3] ;
                for(i in 1:I){
                    ps[1] = log(theta[1]) + von_mises_lpdf(R[i] | ori, kappa1) ;
                    ps[2] = log(theta[2]) + von_mises_lpdf(R[i] | ori, kappa2) ;
                    ps[3] = log(theta[3]) + von_mises_lpdf(R[i] | ori + pi(), kappa3) ;
                    target += D[i] * log_sum_exp(ps) ;
                }
            }

            generated quantities {
                vector[I] log_lik ;

                real ps[3] ;
                for(i in 1:I){
                    ps[1] = log(theta[1]) + von_mises_lpdf(R[i] | ori, kappa1) ;
                    ps[2] = log(theta[2]) + von_mises_lpdf(R[i] | ori, kappa2) ;
                    ps[3] = log(theta[3]) + von_mises_lpdf(R[i] | ori + pi(), kappa3) ;
                    log_lik[i] = D[i] *  log_sum_exp(ps) ;
                }
            }
        """
    elif model == "vonmises":
        model_code = """
            data {
                int I ;
                int D[I] ;
            }

            transformed data {
                real R[I] ;
                for (i in 1:I){
                    R[i] = 2.0 * pi() * i / I ;
                }
            }

            parameters {
                unit_vector[2] O ;
                real<lower=0> kappa ;
            }

            transformed parameters{
                real<lower=-pi(), upper=pi()> ori ;

                // convert unit vector
                ori = atan2(O[1], O[2]) ;
            }

            model {
                for(i in 1:I){
                    target += D[i] * von_mises_lpdf(R[i] | ori, kappa) ;
                }
            }

            generated quantities {
                real MRL ;
                real CV ;
                real CSD ;
                vector[I] log_lik ;

                MRL = modified_bessel_first_kind(1, kappa) / modified_bessel_first_kind(0, kappa) ;
                CV = 1 - MRL ;
                CSD = sqrt(-2 * log(MRL)) ;
                for(i in 1:I){
                    log_lik[i] = D[i] * von_mises_lpdf(R[i] | ori, kappa) ;
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
