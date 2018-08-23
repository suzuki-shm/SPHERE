functions {
    real aevon_mises_lpdf(real theta, real mu, real kappa, real nu){
        return von_mises_lpdf(theta + nu * cos(theta - mu) | mu, kappa) ;
    }

    // Normalizing constant for asymmetry extended von mises distribution
    real aedvon_mises_normalize_constraint(real mu, real kappa, real nu, int N){
        real lp ;

        lp = 0 ;
        for (n in 1:N){
            real theta ;
            theta = -pi() + (2.0 * pi() / N) * n ;
            lp += exp(kappa * cos(theta - mu + nu * cos(theta - mu))) ;
        }
        return log(lp) - log(2 * pi() * modified_bessel_first_kind(0, kappa)) ;
    }

    // log probability density of asymmetry extended von mises distribution
    real aedvon_mises_lpdf(real theta, real mu, real kappa, real nu, int N){
        real logncon ;
        logncon = aedvon_mises_normalize_constraint(mu, kappa, nu, N) ;
        return aevon_mises_lpdf(theta | mu, kappa, nu) - logncon ;
    }

    // log probability density of mixture of asymmetry extended von mises distribution
    real aedvon_mises_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, vector nu, int N) {
        vector[K] lp;
        for (k in 1:K){
            lp[k] = log(a[k]) + aedvon_mises_lpdf(R | mu[k], kappa[k], nu[k], N) ;
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
    vector<lower=0.0>[K] kappa[S] ;
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
        kappa[s] ~ student_t(2.5, 0, 0.2) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * aedvon_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], nu[SUBJECT[i]], L) ;
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
        PTR[s] = exp(2 * kappa[s]) ;
        mPTR[s] = sum(PTR[s] ./ K) ;
        wmPTR[s] = sum(PTR[s] .* alpha) ;
        // Mean resultant length
        for (k in 1:K){
            MRL[s][k] = modified_bessel_first_kind(1, kappa[s][k]) / modified_bessel_first_kind(0, kappa[s][k]) ;
        }
        // Circular variance
        CV[s] = 1.0 - MRL[s] ;
        // Circular standard variation
        CSD[s] = sqrt(-2.0 * log(MRL[s])) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * aedvon_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], nu[SUBJECT[i]], L) ;
    }
}
