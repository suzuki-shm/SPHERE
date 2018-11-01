functions {
    real trans_sin2(real theta, real mu, real nu){
        real theta_mu ;
        theta_mu =  theta - nu * sin(theta - mu) * sin(theta - mu) ;
        return theta_mu ;
    }

    real miaevon_mises_lpdf(real theta, real mu, real kappa, real nu){
        return von_mises_lpdf(trans_sin2(theta, mu, nu)| mu, kappa) ;
    }

    real miaevon_mises_normalize_constraint(real mu, real kappa, real nu, int N){
        // Numerical integration by composite Simpson's rule
        vector[N+1] lp ;
        real h ;

        h = 2 * pi() / N ;
        lp[1] = miaevon_mises_lpdf(-pi() | mu, kappa, nu) ;
        for (n in 1:(N/2)){
            lp[2*n] = log(4) + miaevon_mises_lpdf(-pi() + h*(2*n-1) | mu, kappa, nu) ;
        }
        for (n in 1:(N/2-1)){
            lp[2*n+1] = log(2) + miaevon_mises_lpdf(-pi() + h*2*n | mu, kappa, nu) ;
        }
        lp[N+1] = miaevon_mises_lpdf(pi() | mu, kappa, nu) ;
        return (log(h/3) + log_sum_exp(lp)) ;

    }

    real miaevon_mises_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, vector nu) {
        vector[K] lp;
        real logncon ;

        for (k in 1:K){
            logncon = miaevon_mises_normalize_constraint(mu[k], kappa[k], nu[k], 20) ;
            lp[k] = log(a[k]) + miaevon_mises_lpdf(R | mu[k], kappa[k], nu[k]) - logncon ;
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
}

transformed data {
    real<lower=-pi(), upper=pi()> RADIAN[I] ;
    vector<lower=0.0>[K] A; //hyperparameter for dirichlet distribution

    for (i in 1:I){
        RADIAN[i] = -pi() + (2.0 * pi() / L) * (LOCATION[i] - 1) ;
    }
    for (k in 1:K){
        A[k] = 50 / k ;
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
        kappa[s] ~ student_t(2.5, 0, 0.2025) ;
        nu[s] ~ normal(0, 1) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * miaevon_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], nu[SUBJECT[i]]) ;
    }
}

generated quantities {
    vector<lower=1.0>[K] PTR[S] ;
    vector<lower=1.0>[K] wPTR[S] ;
    vector<lower=1.0>[S] mwPTR ;
    vector<lower=0.0, upper=1.0>[K] MRL[S] ;
    vector<lower=0.0, upper=1.0>[K] CV[S] ;
    vector<lower=0.0>[K] CSD[S] ;
    vector[I] log_lik ;

    for(s in 1:S){
        // Fold change of max p.d.f. to min p.d.f.
        PTR[s] = exp(2 * kappa[s]) ;
        wPTR[s] = exp(2 * alpha .* kappa[s]) ;
        mwPTR[s] = sum(wPTR[s]) ;
        // Mean resultant length
        for (k in 1:K){
            MRL[s][k] = modified_bessel_first_kind(1, kappa[s][k]) / modified_bessel_first_kind(0, kappa[s][k]) ;
        }
        // Circular variance
        CV[s] = 1 - MRL[s] ;
        // Circular standard variation
        CSD[s] = sqrt(-2 * log(MRL[s])) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * miaevon_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], nu[SUBJECT[i]]) ;
    }
}
