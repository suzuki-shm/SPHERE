functions{
    real wrappedcauchy_lpdf(real theta, real mu, real rho){
        return log(1 - pow(rho, 2)) - log(2) - log(pi()) - log(1 + pow(rho, 2) - 2 * rho * cos(theta - mu)) ;
    }

    real dwrappedcauchy_normalize_constraint(real mu, real rho, int N){
        real lp;
        lp = 0 ;
        for (n in 1:N){
            real theta ;
            theta = -pi() + (2.0 * pi() / N) * (n-1) ;
            lp += (1 - pow(rho, 2))  / (2 * pi() * (1 + pow(rho, 2) - 2 * rho * cos(theta - mu))) ;
        }
        return log(lp) ;
    }

    real dwrappedcauchy_lpdf(real theta, real mu, real rho, int N){
        real logncon ;
        logncon = dwrappedcauchy_normalize_constraint(mu, rho, N) ;
        return wrappedcauchy_lpdf(theta | mu, rho) - logncon ;
    }

    real dwrappedcauchy_mixture_lpdf(real R, int K, vector a, vector mu, vector rho, int N){
        vector[K] lp ;
        for (k in 1:K){
            lp[k] = log(a[k]) + dwrappedcauchy_lpdf(R | mu[k], rho[k], N) ;
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
    vector<lower=0.0, upper=1.0>[K] rho[S] ;
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
        rho[s] ~ student_t(2.5, 0, 0.17) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * dwrappedcauchy_mixture_lpdf(RADIAN[i] | K, alpha, ori, rho[SUBJECT[i]], L) ;
    }
}

generated quantities {
    vector<lower=0.0>[K] kappa[S] ;
    vector<lower=1.0>[K] PTR[S] ;
    vector<lower=1.0>[S] mPTR ;
    vector<lower=1.0>[S] wmPTR ;
    vector<lower=0.0, upper=1.0>[K] MRL[S] ;
    vector<lower=0.0, upper=1.0>[K] CV[S] ;
    vector<lower=0.0>[K] CSD[S] ;
    vector[I] log_lik ;

    for(s in 1:S){
        // See (Jones&Pewsey, 2005) about this transformation
        kappa[s] = 2 * atanh(rho[s]) ;
        // Fold change of max p.d.f. to min p.d.f.
        PTR[s] = exp(2 * kappa[s]) ;
        mPTR[s] = sum(PTR[s] ./ K) ;
        wmPTR[s] = sum(PTR[s] .* alpha) ;
        // Mean resultant length
        MRL[s] = rho[s] ;
        // Circular variance
        CV[s] = 1 - MRL[s] ;
        // Circular standard variation
        CSD[s] = sqrt(-2 * log(MRL[s])) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * dwrappedcauchy_mixture_lpdf(RADIAN[i] | K, alpha, ori, rho[SUBJECT[i]], L) ;
    }
}
