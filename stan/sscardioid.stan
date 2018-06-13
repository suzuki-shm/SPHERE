functions{
    real sscardioid_lpdf(real theta, real mu, real rho, real lambda){
        return log(1 + 2 * rho * cos(theta - mu)) - log(2) - log(pi()) + log(1 + lambda * sin(theta - mu)) ;
    }

    real sscardioid_mixture_lpdf(real R, int K, vector a, vector mu, vector rho, vector lambda) {
        vector[K] lp;
        for (k in 1:K){
            lp[k] = log(a[k]) + sscardioid_lpdf(R | mu[k], rho[k], lambda[k]) ;
        }
        return log_sum_exp(lp) ;
    }
}

data {
    int I ;
    int S ;
    int L ;
    int<lower=1, upper=L> LOCATION[I] ;
    int<lower=1, upper=S> SUBJECT[I] ;
    int<lower=0> DEPTH[I] ;
    int<lower=1> K ; // number of mixed distribution
    vector<lower=0>[K] A; //hyperparameter for dirichlet distribution
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
    simplex[K] alpha ;
    unit_vector[2] O[K] ;
    vector<lower=0, upper=0.5>[K] rho[S] ;
    vector<lower=-1.0, upper=1.0>[K] lambda[S] ;
}

transformed parameters{
    vector[K] ori ;

    // convert unit vector
    for (k in 1:K){
        ori[k] = atan2(O[k][1], O[k][2]) ;
    }
}

model {
    alpha ~ dirichlet(A) ;
    for(s in 1:S){
        rho[s] ~ normal(0.25, 0.25) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * sscardioid_mixture_lpdf(RADIAN[i] | L, alpha, ori, rho[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
}

generated quantities {
    vector[K] PTR[S] ;
    vector[K] MRL[S] ;
    vector[K] CV[S] ;
    vector[K] CSD[S] ;
    vector[I] log_lik ;

    for(s in 1:S){
        // Fold change of max p.d.f. to min p.d.f.
        PTR[s] = (1 + 2 * rho[s]) ./ (1 - 2 * rho[s]) ;
        // Mean resultant length
        MRL[s] = rho[s] ;
        // Circular variance
        CV[s] = 1 - MRL[s] ;
        // Circular standard variation
        CSD[s] = sqrt(-2 * log(MRL[s])) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * sscardioid_mixture_lpdf(RADIAN[i] | K, alpha, ori, rho[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
}
