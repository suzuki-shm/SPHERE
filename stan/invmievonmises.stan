functions {
    real inv_trans_APF(real theta, real mu, real nu, real lambda){
        real t ;
        t = theta ;
        for (i in 1:8){
            t = t - ((t-mu - nu * sin(t-mu) + lambda * pow(sin(t-mu - nu * sin(t-mu)), 2) - theta) / ((1 + 2 * lambda * sin(t-mu - nu  * sin(t-mu)) * cos(t-mu - nu * sin(t-mu)))* (1 - nu * cos(t-mu)))) ;
        }
        return t ;
    }

    real invmievon_mises_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, vector nu, vector lambda) {
        vector[K] lp;
        real alpha1 ;
        for (k in 1:K){
            alpha1 = modified_bessel_first_kind(1, kappa[k]) / modified_bessel_first_kind(0, kappa[k]) ;
            lp[k] = log(a[k]) + von_mises_lpdf(inv_trans_APF(R, mu[k], nu[k], lambda[k]) | mu[k], kappa[k]) - log(1 - nu[k] * alpha1) ;
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
    // peakness parameter as raw
    vector<lower=-1.0, upper=1.0>[K] lambda_raw[S] ;
}

transformed parameters{
    vector[K] ori ;
    // peakness parameter
    vector<lower=-1.0, upper=1.0>[K] lambda[S] ;

    // convert unit vector
    for (k in 1:K){
        ori[k] = atan2(O[k][1], O[k][2]) ;
    }
    for (s in 1:S){
        lambda[s] = lambda_raw[s] .* (1 - fabs(nu[s])) ;
    }
}

model {
    alpha ~ dirichlet(A) ;
    for(s in 1:S){
        kappa[s] ~ student_t(2.5, 0, 0.2) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * invmievon_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], nu[SUBJECT[i]], lambda[SUBJECT[i]]) ;
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
        log_lik[i] = DEPTH[i] * invmievon_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], nu[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
}
