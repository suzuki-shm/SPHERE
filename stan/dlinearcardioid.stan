 functions {
    real linearcardioid_lpdf(real theta, real mu, real rho){
        return log(1 + 2 * rho * (fabs(fabs(theta - mu) - pi()) - pi() / 2)) - log(2) -log(pi())   ;
    }

    real dlinearcardioid_normalize_constraint(int N){
        return (N + 2) / (2 * pi()) ;
    }

    real dlinearcardioid_lpdf(real theta, real mu, real rho, int N){
        real logncon ;
        logncon = dlinearcardioid_normalize_constraint(N) ;
        return linearcardioid_lpdf(theta | mu, rho) - logncon ;
    }

    real dlinearcardioid_mixture_lpdf(real R, int K, vector a, vector mu, vector rho, int N) {
        vector[K] lp;
        for (k in 1:K){
            lp[k] = log(a[k]) + dlinearcardioid_lpdf(R | mu[k], rho[k], N) ;
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
    vector<lower=0.0, upper=1/pi()>[K] rho[S] ;
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
        rho[s] ~ student_t(2.5, 0, 0.105) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * dlinearcardioid_mixture_lpdf(RADIAN[i]| K, alpha, ori, rho[SUBJECT[i]], L) ;
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
        PTR[s] = exp(2 * rho[s]) ;
        mPTR[s] = sum(PTR[s] ./ K) ;
        wmPTR[s] = sum(PTR[s] .* alpha) ;
        MRL[s] = rho[s] ;
        CV[s] = 1 - MRL[s] ;
        CSD[s] = sqrt(-2 * log(rho[s])) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * dlinearcardioid_mixture_lpdf(RADIAN[i]| K, alpha, ori, rho[SUBJECT[i]], L) ;
    }
}
