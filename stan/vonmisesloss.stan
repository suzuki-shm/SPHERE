functions {
    real von_mises_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, int L) {
        vector[K] lp;
        for (k in 1:K){
            lp[k] = log(a[k]) + von_mises_lpdf(R | mu[k], kappa[k]) + log(2.0) + log(pi()) - log(L) ;
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
    vector<lower=-pi(), upper=pi()>[I] RADIAN ;
    vector<lower=0.0>[K] A; //hyperparameter for dirichlet distribution
    int DEPTH_SUM[S] ;

    RADIAN = -pi() + (2.0 * pi() / L) * (to_vector(LOCATION) - 1) ;
    A = rep_vector(50.0/K, K) ;
    for (s in 1:S){
        DEPTH_SUM[s] = 0 ;
    }
    for (i in 1:I){
        DEPTH_SUM[SUBJECT[i]] += DEPTH[i] ;
    }
}

parameters {
    simplex[K] alpha ;
    unit_vector[2] O[K] ;
    // Unconstrained concentration parameter
    vector[K] kappa_uncon[S] ;
}

transformed parameters{
    vector[K] ori ;
    vector<lower=0>[K] kappa[S] ;

    // convert unit vector
    for (k in 1:K){
        ori[k] = atan2(O[k][1], O[k][2]) ;
    }
    // Add upper bound to kappa using alpha (see 'Lower and Upper Bounded Scalar' in Stan manual)
    for (s in 1:S){
        kappa[s] = log(4) ./ (2 * alpha) .* inv_logit(kappa_uncon[s]) ;
    }
}

model {
    alpha ~ dirichlet(A) ;
    for (s in 1:S){
        alpha .* kappa[s] ~ student_t(2.5, 0, 0.2025) ;
        // Jacobian adjustment for parameter transformation (see 'Lower and Upper Bounded Scalar' in Stan manual)
        target += log(log(4) ./ (2 * alpha)) + log_inv_logit(kappa_uncon[s]) + log1m_inv_logit(kappa_uncon[s]) ;
    }
    for(i in 1:I){
        DEPTH[i] ~ poisson(DEPTH_SUM[SUBJECT[i]] * exp(von_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], L))) ;
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
    real log_lik_sum ;

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
        log_lik[i] = poisson_lpmf(DEPTH[i] | DEPTH_SUM[SUBJECT[i]] * exp(von_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], L))) ;
    }
    log_lik_sum = sum(log_lik) ;
}
