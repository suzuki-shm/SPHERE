 functions {
    real linearcardioid_lpdf(real theta, real mu, real rho){
        return log(1 + 2 * rho * (fabs(fabs(theta - mu) - pi()) - pi() / 2)) - log(2) -log(pi())   ;
    }

    real linearcardioid_mixture_lpdf(real R, int K, vector a, vector mu, vector rho, int L) {
        vector[K] lp;
        for (k in 1:K){
            lp[k] = log(a[k]) + linearcardioid_lpdf(R | mu[k], rho[k]) + log(2.0) + log(pi()) - log(L) ;
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
    vector[K] rho_uncon[S] ;
}

transformed parameters{
    vector[K] ori ;
    vector<lower=0.0, upper=1.0/pi()>[K] rho[S] ;

    // convert unit vector
    for (k in 1:K){
        ori[k] = atan2(O[k][1], O[k][2]) ;
    }
    // Add upper bound to kappa using alpha (see 'Lower and Upper Bounded Scalar' in Stan manual)
    for (s in 1:S){
        for (k in 1:K){
            rho[s][k] = fmin(3.0/5.0/pi()/alpha[k], 1.0/pi()) .* inv_logit(rho_uncon[s][k]) ;
        }
    }
}

model {
    alpha ~ dirichlet(A) ;
    for(s in 1:S){
        alpha .* rho[s] ~ student_t(2.5, 0, 0.1) ;
        // Jacobian adjustment for parameter transformation (see 'Lower and Upper Bounded Scalar' in Stan manual)
        for (k in 1:K){
            target += log(fmin(3.0/5.0/pi()/alpha[k], 1.0/pi())) + log_inv_logit(rho_uncon[s][k]) + log1m_inv_logit(rho_uncon[s][k]) ;
        }
        // Jacobian adjustment for alpha * concentration parameter
        target += -log(alpha) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * linearcardioid_mixture_lpdf(RADIAN[i]| K, alpha, ori, rho[SUBJECT[i]], L) ;
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
        PTR[s] = (1 + pi() * rho[s])  ./ (1 - pi() * rho[s]) ;
        wPTR[s] = (1 + pi() * alpha .* rho[s])  ./ (1 - pi() * alpha .* rho[s]) ;
        mwPTR[s] = mean(wPTR[s]) ;
        MRL[s] = rho[s] ;
        CV[s] = 1 - MRL[s] ;
        CSD[s] = sqrt(-2 * log(rho[s])) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * linearcardioid_mixture_lpdf(RADIAN[i]| K, alpha, ori, rho[SUBJECT[i]], L) ;
    }
    log_lik_sum = sum(log_lik) ;
}
