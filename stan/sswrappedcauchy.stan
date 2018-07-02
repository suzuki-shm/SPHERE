functions{
    real sswrappedcauchy_lpdf(real theta, real mu, real kappa, real lambda){
        return log(1 - pow(kappa, 2)) - log(2) - log(pi()) - log(1 + pow(kappa, 2) - 2 * kappa * cos(theta - mu)) + log(1 + lambda * sin(theta - mu)) ;
    }

    real sswrappedcauchy_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, vector lambda) {
        vector[K] lp;
        for (k in 1:K){
            lp[k] = log(a[k]) + sswrappedcauchy_lpdf(R | mu[k], kappa[k], lambda[k]) ;
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
    vector<lower=0, upper=1.0>[K] kappa[S] ;
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
        kappa[s] ~ normal(0.5, 0.5) ;
        lambda[s] ~ normal(0, 1) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * sswrappedcauchy_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
}

generated quantities {
    vector<lower=1.0>[K] PTR[S] ;
    vector<lower=1.0>[S] mPTR ;
    vector<lower=1.0>[S] wPTR ;
    vector<lower=0.0, upper=1.0>[K] MRL[S] ;
    vector<lower=0.0, upper=1.0>[K] CV[S] ;
    vector<lower=0.0>[K] CSD[S] ;
    vector[I] log_lik ;

    for(s in 1:S){
        // Fold change of max p.d.f. to min p.d.f.
        PTR[s] = (1 + kappa[s] .* kappa[s]) ./ ((1 - kappa[s]) .* (1 - kappa[s])) ;
        mPTR[s] = mean((1 + kappa[s] .* kappa[s] / K) ./ ((1 - kappa[s] / K) .* (1 - kappa[s] / K))) ;
        wPTR[s] = mean((1 + kappa[s] .* kappa[s] .* alpha) ./ ((1 - kappa[s] .* alpha) .* (1 - kappa[s] .* alpha))) ;
        // Mean resultant length
        MRL[s] = kappa[s] ;
        // Circular variance
        CV[s] = 1 - MRL[s] ;
        // Circular standard variation
        CSD[s] = sqrt(-2 * log(MRL[s])) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * sswrappedcauchy_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
}
