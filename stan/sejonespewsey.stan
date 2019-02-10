functions {
    real jonespewsey_lpdf(real theta, real mu, real kappa, real psi){
        if (fabs(psi) < 1e-10){
            return kappa * cos(theta - mu) ;
        }else{
            return 1 / psi * log(cosh(kappa * psi) + sinh(kappa * psi)*cos(theta - mu)) ;
        }
    }

    real sevon_mises_lpdf(real theta, real mu, real kappa, real lambda){
        return von_mises_lpdf(theta + lambda * sin(theta-mu) | mu, kappa) ;
    }

    real sejonespewsey_lpdf(real theta, real mu, real kappa, real psi, real lambda){
        return jonespewsey_lpdf(theta + lambda * sin(theta-mu) | mu, kappa, psi) ;
    }

    real sejonespewsey_normalize_constraint(real mu, real kappa, real psi, real lambda, int N){
        // Numerical integration by composite Simpson's rule
        vector[N+1] lp ;
        real h ;
        real logncon ;
        h = 2 * pi() / N ;

        if (fabs(psi) < 1e-10){
            lp[1] = sevon_mises_lpdf(-pi() | mu, kappa, lambda) ;
            for (n in 1:(N/2)){
                lp[2*n] = log(4) + sevon_mises_lpdf(-pi() + h*(2*n-1) | mu, kappa, lambda) ;
            }
            for (n in 1:(N/2-1)){
                lp[2*n+1] = log(2) + sevon_mises_lpdf(-pi() + h*2*n | mu, kappa, lambda) ;
            }
            lp[N+1] = sevon_mises_lpdf(pi() | mu, kappa, lambda) ;
        }else{
            lp[1] = sejonespewsey_lpdf(-pi() | mu, kappa, psi, lambda) ;
            for (n in 1:(N/2)){
                lp[2*n] = log(4) + sejonespewsey_lpdf(-pi() + h*(2*n-1) | mu, kappa, psi, lambda) ;
            }
            for (n in 1:(N/2-1)){
                lp[2*n+1] = log(2) + sejonespewsey_lpdf(-pi() + h*2*n | mu, kappa, psi, lambda) ;
            }
            lp[N+1] = sejonespewsey_lpdf(pi() | mu, kappa, psi, lambda) ;
        }
        logncon =  log(h/3) + log_sum_exp(lp) ;
        return logncon ;
    }

    real sejonespewsey_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, vector psi, vector lambda){
        vector[K] lp ;
        real logncon ;

        for (k in 1:K){
            logncon = sejonespewsey_normalize_constraint(mu[k], kappa[k], psi[k], lambda[k], 20) ;
            lp[k] = log(a[k]) + sejonespewsey_lpdf(R | mu[k], kappa[k], psi[k], lambda[k]) - logncon ;
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
    vector<lower=-1.0, upper=1.0>[K] psi[S] ;
    // peakness parameter
    vector<lower=-1.0, upper=1.0>[K] lambda[S] ;
}

transformed parameters{
    vector[K] ori ;
    vector<lower=0.0>[K] kappa[S] ;

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
        psi[s] ~ normal(0, 1) ;
        lambda[s] ~ normal(0, 1) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * sejonespewsey_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], psi[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
}

generated quantities {
    vector<lower=1.0>[K] PTR[S] ;
    vector<lower=1.0>[K] wPTR[S] ;
    vector<lower=1.0>[S] mwPTR ;
    vector[I] log_lik ;

    for(s in 1:S){
        // Fold change of max p.d.f. to min p.d.f.
        PTR[s] = exp(2 * kappa[s]) ;
        wPTR[s] = exp(2 * alpha .* kappa[s]) ;
        mwPTR[s] = sum(wPTR[s]) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * sejonespewsey_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], psi[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
}
