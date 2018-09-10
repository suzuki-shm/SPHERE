functions {
    real inv_trans_sin2(real theta, real mu, real nu){
        real t ;
        t = theta ;
        for (i in 1:8){
            t = t - ((t + nu*pow(sin(t-mu),2) - theta) / (1 + 2 * nu * sin(t-mu) * cos(t-mu))) ;
        }
        return t ;
    }

    real jonespewsey_lpdf(real theta, real mu, real kappa, real psi){
        if (fabs(psi) < 1e-10){
            return kappa * cos(theta - mu) ;
        }else{
            return 1 / psi * log(cosh(kappa * psi) + sinh(kappa * psi)*cos(theta - mu)) ;
        }
    }

    real jonespewsey_normalize_constraint(real mu, real kappa, real psi, int N){
        // Numerical integration by composite Simpson's rule
        vector[N+1] lp ;
        real h ;
        real logncon ;

        if (fabs(psi) < 1e-10){
            logncon = log(2 * pi() * modified_bessel_first_kind(0, kappa)) ;
        }else{
            h = 2 * pi() / N ;
            lp[1] = jonespewsey_lpdf(-pi() | mu, kappa, psi) ;
            for (n in 1:(N/2)){
                lp[2*n] = log(4) + jonespewsey_lpdf(-pi() + h*(2*n-1) | mu, kappa, psi) ;
            }
            for (n in 1:(N/2-1)){
                lp[2*n+1] = log(2) + jonespewsey_lpdf(-pi() + h*2*n | mu, kappa, psi) ;
            }
            lp[N+1] = jonespewsey_lpdf(pi() | mu, kappa, psi) ;
            logncon =  log(h/3) + log_sum_exp(lp) ;
        }
        return logncon ;
    }

    real invmijonespewsey_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, vector psi, vector nu){
        vector[K] lp ;
        real logncon ;

        for (k in 1:K){
            logncon = jonespewsey_normalize_constraint(mu[k], kappa[k], psi[k], 20) ;
            lp[k] = log(a[k]) + jonespewsey_lpdf(inv_trans_sin2(R, mu[k], nu[k]) | mu[k], kappa[k], psi[k]) - logncon ;
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
    vector<lower=-1.0, upper=1.0>[K] psi[S] ;
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
        target += DEPTH[i] * invmijonespewsey_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], psi[SUBJECT[i]], nu[SUBJECT[i]]) ;
    }
}

generated quantities {
    vector<lower=1.0>[K] PTR[S] ;
    vector<lower=1.0>[S] mPTR ;
    vector<lower=1.0>[S] wmPTR ;
    vector[I] log_lik ;

    for(s in 1:S){
        // Fold change of max p.d.f. to min p.d.f.
        PTR[s] = exp(2 * kappa[s]) ;
        mPTR[s] = sum(PTR[s] ./ K) ;
        wmPTR[s] = sum(PTR[s] .* alpha) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * invmijonespewsey_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], psi[SUBJECT[i]], nu[SUBJECT[i]]) ;
    }
}
