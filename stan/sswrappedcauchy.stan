functions{
    real sswrappedcauchy_lpdf(real theta, real mu, real rho, real lambda){
        return log(1 - pow(rho, 2)) - log(2) - log(pi()) - log(1 + pow(rho, 2) - 2 * rho * cos(theta - mu)) + log(1 + lambda * sin(theta - mu)) ;
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
        if(i < L/2.0){
            RADIAN[i] = 2.0 * pi() * LOCATION[i] / L ;
        }else{
            RADIAN[i] = 2.0 * pi() * (LOCATION[i] - L) / L ;
        }
    }
}

parameters {
    unit_vector[2] O ;
    real<lower=0, upper=1.0> rho[S] ;
    real<lower=-1.0, upper=1.0> lambda[S] ;
}

transformed parameters{
    real<lower=-pi(), upper=pi()> ori ;

    // convert unit vector
    ori = atan2(O[1], O[2]) ;
}

model {
    for(i in 1:I){
        target += DEPTH[i] * sswrappedcauchy_lpdf(RADIAN[i] | ori, rho[SUBJECT[i]], lambda[SUBJECT[i]]) ;
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
        log_lik[i] = DEPTH[i] * sswrappedcauchy_lpdf(RADIAN[i] | ori, rho[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
}
