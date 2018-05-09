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
        if(i < L/2){
            RADIAN[i] = 2.0 * pi() * LOCATION[i] / L ;
        }else{
            RADIAN[i] = 2.0 * pi() * (LOCATION[i] - L) / L ;
        }
    }
}

parameters {
    unit_vector[2] O ;
    real<lower=0, upper=1/pi()> rho[S] ;
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
