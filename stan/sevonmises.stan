functions {
    real sevon_mises_lpdf(real theta, real mu, real kappa, real lambda){
        return von_mises_lpdf(theta + lambda * sin(theta-mu) | mu, kappa) ;
    }

    real sevon_mises_normalize_constraint(real mu, real kappa, real lambda, int N){
        // numerical integration by composite Simpson's rule
        vector[N+1] lp ;
        real h ;

        h = 2 * pi() / N ;
        lp[1] = sevon_mises_lpdf(-pi() | mu, kappa, lambda) ;
        for (n in 1:(N/2)){
            lp[2*n] = log(4) + sevon_mises_lpdf(-pi() + h*(2*n-1) | mu, kappa, lambda) ;
        }
        for (n in 1:(N/2-1)){
            lp[2*n+1] = log(2) + sevon_mises_lpdf(-pi() + h*2*n | mu, kappa, lambda) ;
        }
        lp[N+1] = sevon_mises_lpdf(pi() | mu, kappa, lambda) ;
        return (log(h/3) + log_sum_exp(lp)) ;

    }

    real sevon_mises_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, vector lambda) {
        vector[K] lp;
        real logncon ;

        for (k in 1:K){
            logncon = sevon_mises_normalize_constraint(mu[k], kappa[k], lambda[k], 20) ;
            lp[k] = log(a[k]) + sevon_mises_lpdf(R | mu[k], kappa[k], lambda[k]) - logncon ;
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

    RADIAN = -pi() + (2.0 * pi() / L) * (to_vector(LOCATION) - 1) ;
    A = rep_vector(50.0/K, K) ;
}

parameters {
    simplex[K] alpha ;
    unit_vector[2] O[K] ;
    // Unconstrained concentration parameter
    vector[K] kappa_uncon[S] ;
    // peakness parameter
    vector<lower=-1.0, upper=1.0>[K] lambda[S] ;
}

transformed parameters{
    vector[K] ori ;
    vector<lower=0, upper=4.0>[K] kappa[S] ;

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
    for(s in 1:S){
        alpha .* kappa[s] ~ student_t(2.5, 0, 0.2025) ;
        // Jacobian adjustment for parameter transformation (see 'Lower and Upper Bounded Scalar' in Stan manual)
        target += log(log(4) ./ (2 * alpha)) + log_inv_logit(kappa_uncon[s]) + log1m_inv_logit(kappa_uncon[s]) ;
        lambda[s] ~ normal(0, 1) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * sevon_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
}

generated quantities {
    vector<lower=1.0>[K] PTR[S] ;
    vector<lower=1.0>[K] wPTR[S] ;
    vector<lower=1.0>[S] mwPTR ;
    vector[I] log_lik ;
    real log_lik_sum ;

    for(s in 1:S){
        // Fold change of max p.d.f. to min p.d.f.
        PTR[s] = exp(2 * kappa[s]) ;
        wPTR[s] = exp(2 * alpha .* kappa[s]) ;
        mwPTR[s] = mean(wPTR[s]) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * sevon_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
    log_lik_sum = sum(log_lik) ;
}
