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
        if(i < L/2){
            RADIAN[i] = 2.0 * pi() * LOCATION[i] / L ;
        }else{
            RADIAN[i] = 2.0 * pi() * (LOCATION[i] - L) / L ;
        }
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
