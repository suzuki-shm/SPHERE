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
                        choices=["linearcardioid",
                                 "cardioid",
                                 "wrappedcauchy",
                                 "vonmises"],
                        default="vonmises",
                        help="file path of compiled stan model")
    parser.set_defaults(trans=False)
    args = parser.parse_args()
    return vars(args)


def compile_model(output_path=None, model="vonmises"):
    # Stanのモデルを読み込んでコンパイルする
    if model == "linearcardioid":
        model_code = """
            functions {
                real linearcardioid_lpdf(real theta, real mu, real rho){
                    return log(1 / (2 * pi()) * (1 + 2 * rho * (fabs(fabs(theta - mu) - pi()) - pi() / 2))) ;
                }
            }
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
                real<lower=0, upper=0.5> rho[S] ;
                real<lower=0> sigma_rho ;
            }

            transformed parameters{
                real<lower=-pi(), upper=pi()> ori ;
                // convert unit vector
                ori = atan2(O[1], O[2]) ;
            }

            model {
                for(s in 1:S){
                    rho[s] ~ normal(0, sigma_rho) ;
                }
                for(i in 1:I){
                    target += DEPTH[i] * linearcardioid_lpdf(RADIAN[i]| ori, rho[SUBJECT[i]]) ;
                }
            }

            generated quantities {
                vector[I] log_lik ;
                real<lower=0.0, upper=1.0> MRL[S] ;
                real<lower=0.0, upper=1.0> CV[S] ;
                real<lower=0> CSD[S] ;
                real<lower=1.0> PTR[S] ;

                for(s in 1:S){
                    PTR[s] = (1 + pi() * rho[s]) / (1 - pi() * rho[s]) ;
                    MRL[s] = rho[s] ;
                    CV[s] = 1 - MRL[s] ;
                    CSD[s] = sqrt(-2 * log(rho[s])) ;
                }
                for(i in 1:I){
                    log_lik[i] = DEPTH[i] * linearcardioid_lpdf(RADIAN[i]| ori, rho[SUBJECT[i]]) ;
                }
            }
        """
    elif model == "cardioid":
        model_code = """
            functions{
                real cardioid_lpdf(real theta, real mu, real rho){
                    return log(1 / (2 * pi()) * (1 + 2 * rho * cos(theta - mu))) ;
                }
            }
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
                real<lower=0, upper=0.5> rho[S] ;
                real<lower=0> sigma_rho;
            }

            transformed parameters{
                real<lower=-pi(), upper=pi()> ori ;

                // convert unit vector
                ori = atan2(O[1], O[2]) ;
            }

            model {
                for(s in 1:S){
                    rho[s] ~ normal(0, sigma_rho) ;
                }
                for(i in 1:I){
                    target += DEPTH[i] * cardioid_lpdf(RADIAN[i] | ori, rho[SUBJECT[i]]) ;
                }
            }

            generated quantities {
                real<lower=1.0> PTR[S] ;
                real MRL[S] ;
                real CV[S] ;
                real CSD[S] ;
                vector[I] log_lik ;

                for(s in 1:S){
                    // Fold change of max p.d.f. to min p.d.f.
                    PTR[s] = (1 + 2 * rho[s]) / (1 - 2 * rho[s]) ;
                    // Mean resultant length
                    MRL[s] = rho[s] ;
                    // Circular variance
                    CV[s] = 1 - MRL[s] ;
                    // Circular standard variation
                    CSD[s] = sqrt(-2 * log(MRL[s])) ;
                }
                for(i in 1:I){
                    log_lik[i] = DEPTH[i] * cardioid_lpdf(RADIAN[i] | ori, rho[SUBJECT[i]]) ;
                }
            }
        """
    elif model == "wrappedcauchy":
        model_code = """
            functions{
                real wrappedcauchy_lpdf(real theta, real mu, real rho){
                    return log((1 - pow(rho, 2)) / (2 * pi() * (1 + pow(rho, 2) - 2 * rho * cos(theta - mu)))) ;
                }
            }
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
                real<lower=0, upper=1.0> rho[S] ;
                real<lower=0> sigma_rho;
            }

            transformed parameters{
                real<lower=-pi(), upper=pi()> ori ;

                // convert unit vector
                ori = atan2(O[1], O[2]) ;
            }

            model {
                for(s in 1:S){
                    rho[s] ~ normal(0, sigma_rho) ;
                }
                for(i in 1:I){
                    target += DEPTH[i] * wrappedcauchy_lpdf(RADIAN[i] | ori, rho[SUBJECT[i]]) ;
                }
            }

            generated quantities {
                real<lower=1.0> PTR[S] ;
                real MRL[S] ;
                real CV[S] ;
                real CSD[S] ;
                vector[I] log_lik ;

                for(s in 1:S){
                    // Fold change of max p.d.f. to min p.d.f.
                    PTR[s] = (1 + pow(rho[s], 2)) / pow(1 - rho[s], 2) ;
                    // Mean resultant length
                    MRL[s] = rho[s] ;
                    // Circular variance
                    CV[s] = 1 - MRL[s] ;
                    // Circular standard variation
                    CSD[s] = sqrt(-2 * log(MRL[s])) ;
                }
                for(i in 1:I){
                    log_lik[i] = DEPTH[i] * wrappedcauchy_lpdf(RADIAN[i] | ori, rho[SUBJECT[i]]) ;
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
                real<lower=0> sigma_kappa ;
            }

            transformed parameters{
                real<lower=-pi(), upper=pi()> ori ;

                // convert unit vector
                ori = atan2(O[1], O[2]) ;
            }

            model {
                for(s in 1:S){
                    kappa[s] ~ normal(0, sigma_kappa) ;
                }
                for(i in 1:I){
                    target += DEPTH[i] * von_mises_lpdf(RADIAN[i] | ori, kappa[SUBJECT[i]]) ;
                }
            }

            generated quantities {
                real<lower=1.0> PTR[S] ;
                real MRL[S] ;
                real CV[S] ;
                real CSD[S] ;
                vector[I] log_lik ;

                for(s in 1:S){
                    // Fold change of max p.d.f. to min p.d.f.
                    PTR[s] = exp(2 * kappa[s]) ;
                    // Mean resultant length
                    MRL[s] = modified_bessel_first_kind(1, kappa[s]) / modified_bessel_first_kind(0, kappa[s]) ;
                    // Circular variance
                    CV[s] = 1 - MRL[s] ;
                    // Circular standard variation
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
