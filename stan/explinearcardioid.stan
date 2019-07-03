 functions {
    real explinearcardioid_lpdf(real theta, real mu, real rho){
        return 2 * rho * (fabs(fabs(theta - mu) - pi()) - pi() / 2) + log(rho) - log(exp(pi()*rho) - exp(-pi()*rho))   ;
    }

    real explinearcardioid_mixture_lpdf(real R, int K, vector a, vector mu, vector rho) {
        vector[K] lp;
        for (k in 1:K){
            lp[k] = log(a[k]) + explinearcardioid_lpdf(R | mu[k], rho[k]) ;
        }
        return log_sum_exp(lp) ;
    }

    real uniform_loglik_Simpson(real lower, real upper){
        return log((upper - lower) / 2 / pi()) ;
    }
    
    real explinearcardioid_loglik_Simpson(real lower, real upper, int K, vector a, vector mu, vector rho){
        int M = 20;
        vector[M+1] lp;
        real h;
        h = (upper - lower) / M ;
        lp[1] = explinearcardioid_mixture_lpdf(lower | K, a, mu, rho) ;
        for (m in 1:M/2){
            lp[2*m] = log(4) + explinearcardioid_mixture_lpdf(lower + h*(2*m-1) | K, a, mu, rho) ;
        }
        for (m in 1:M/2-1){
            lp[2*m+1] = log(2) + explinearcardioid_mixture_lpdf(lower + h*2*m | K, a, mu, rho) ;
        }
        lp[M+1] = explinearcardioid_mixture_lpdf(upper | K, a, mu, rho) ;
        return (log(h/3) + log_sum_exp(lp)) ;
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
    real S_uniform ;

    RADIAN = -pi() + (2.0 * pi() / L) * (to_vector(LOCATION) - 1) ;
    A = rep_vector(50.0/K, K) ;
    S_uniform = uniform_loglik_Simpson(-pi(), pi()) ;
}

parameters {
    simplex[K] alpha ;
    unit_vector[2] O[K] ;
    // Unconstrained concentration parameter
    vector[K] rho_uncon[S] ;
}

transformed parameters{
    vector[K] ori ;
    vector<lower=0.0>[K] rho[S] ;
    real S_explinearcardioid[S] ;

    // convert unit vector
    for (k in 1:K){
        ori[k] = atan2(O[k][1], O[k][2]) ;
    }
    // Add upper bound to kappa using alpha (see 'Lower and Upper Bounded Scalar' in Stan manual)
    for (s in 1:S){
        for (k in 1:K){
            rho[s][k] = 1.0/2.0/pi()/alpha[k]*log(4) .* inv_logit(rho_uncon[s][k]) ;
        }
        S_explinearcardioid[s] = explinearcardioid_loglik_Simpson(-pi(), pi(), K, alpha, ori, rho[s]) ;
    }
}

model {
    alpha ~ dirichlet(A) ;
    for(s in 1:S){
        alpha .* rho[s] ~ student_t(2.5, 0, 0.1103) ;
        // Jacobian adjustment for parameter transformation (see 'Lower and Upper Bounded Scalar' in Stan manual)
        for (k in 1:K){
            target += 1.0/2.0/pi()/alpha[k]*log(4) + log_inv_logit(rho_uncon[s][k]) + log1m_inv_logit(rho_uncon[s][k]) ;
        }
    }
    for(i in 1:I){
        target += DEPTH[i] * (explinearcardioid_mixture_lpdf(RADIAN[i]| K, alpha, ori, rho[SUBJECT[i]]) + log(2 * pi()) - log(L) + S_uniform - S_explinearcardioid[SUBJECT[i]]) ;
    }
}

generated quantities {
    vector<lower=1.0>[K] PTR[S] ;
    vector<lower=1.0>[K] wPTR[S] ;
    vector<lower=1.0>[S] mwPTR ;
    vector[I] log_lik ;
    real log_lik_sum ;

    for(s in 1:S){
        PTR[s] = exp(2.0 * pi() * rho[s]) ;
        wPTR[s] = exp(2.0 * pi() * alpha .* rho[s]) ;
        mwPTR[s] = sum(wPTR[s]) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * (explinearcardioid_mixture_lpdf(RADIAN[i]| K, alpha, ori, rho[SUBJECT[i]]) + log(2 * pi()) - log(L) + S_uniform - S_explinearcardioid[SUBJECT[i]]) ;
    }
    log_lik_sum = sum(log_lik) ;
}
