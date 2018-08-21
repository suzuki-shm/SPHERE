functions{
    real wrappedcauchy_lpdf(real theta, real mu, real kappa){
        return log(1 - pow(kappa, 2)) - log(2) - log(pi()) - log(1 + pow(kappa, 2) - 2 * kappa * cos(theta - mu)) ;
    }

    real aewrappedcauchy_lpdf(real theta, real mu, real kappa, real nu){
        return wrappedcauchy_lpdf(theta + nu * cos(theta - mu) | mu, kappa) ;
    }

    real aedwrappedcauchy_normalize_constraint(real mu, real kappa, real nu, int N){
        vector[N] lp ;
        for (n in 1:N){
            real theta ;
            theta = -pi() + (2.0 * pi() / N) * n ;
            lp[n] = aewrappedcauchy_lpdf(theta | mu, kappa, nu) ;
        }
        return log_sum_exp(lp) ;
    }

    real aedwrappedcauchy_lpdf(real theta, real mu, real kappa, real nu, int N){
        real logncon ;
        logncon = aedwrappedcauchy_normalize_constraint(mu, kappa, nu, N) ;
        return aewrappedcauchy_lpdf(theta | mu, kappa, nu) - logncon ;
    }

    real aedwrappedcauchy_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, vector nu, int N){
        vector[K] lp ;
        for (k in 1:K){
            lp[k] = log(a[k]) + aedwrappedcauchy_lpdf(R | mu[k], kappa[k], nu[k], N) ;
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
    vector<lower=0.0>[K] A; //hyperparameter for dirichlet distribution
}

transformed data {
    real<lower=-pi(), upper=pi()> RADIAN[I] ;

    for (i in 1:I){
        RADIAN[i] = -pi() + (2.0 * pi() / L) * (LOCATION[i] - 1) ;
    }
}

parameters {
    simplex[K] alpha ;
    unit_vector[2] O[K] ;
    vector<lower=0.0, upper=1.0>[K] kappa[S] ;
    // skewness parameter
    vector<lower=-1.0, upper=1.0>[K] nu[S] ;
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
        kappa[s] ~ student_t(2.5, 0, 0.17) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * aedwrappedcauchy_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], nu[SUBJECT[i]], L) ;
    }
}

generated quantities {
    vector<lower=1.0>[K] PTR[S] ;
    vector<lower=1.0>[S] mPTR ;
    vector<lower=1.0>[S] wmPTR ;
    vector<lower=0.0, upper=1.0>[K] MRL[S] ;
    vector<lower=0.0, upper=1.0>[K] CV[S] ;
    vector<lower=0.0>[K] CSD[S] ;
    vector[I] log_lik ;

    for(s in 1:S){
        // Fold change of max p.d.f. to min p.d.f.
        PTR[s] = (1 + kappa[s] .* kappa[s]) ./ ((1 - kappa[s]) .* (1 - kappa[s])) ;
        mPTR[s] = sum(PTR[s] ./ K) ;
        wmPTR[s] = sum(PTR[s] .* alpha) ;
        // Mean resultant length
        MRL[s] = kappa[s] ;
        // Circular variance
        CV[s] = 1 - MRL[s] ;
        // Circular standard variation
        CSD[s] = sqrt(-2 * log(MRL[s])) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * aedwrappedcauchy_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], nu[SUBJECT[i]], L) ;
    }
}
