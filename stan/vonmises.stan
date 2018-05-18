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
        if(i < L/2.0){
            RADIAN[i] = 2.0 * pi() * LOCATION[i] / L ;
        }else{
            RADIAN[i] = 2.0 * pi() * (LOCATION[i] - L) / L ;
        }
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
