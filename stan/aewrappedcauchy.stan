functions{
    real wrappedcauchy_lpdf(real theta, real mu, real rho){
        return log(1 - pow(rho, 2)) - log(2) - log(pi()) - log(1 + pow(rho, 2) - 2 * rho * cos(theta - mu)) ;
    }

    real aewrappedcauchy_lpdf(real theta, real mu, real rho, real nu){
        return wrappedcauchy_lpdf(theta + nu * cos(theta - mu) | mu, rho) ;
    }

    real aewrappedcauchy_normalize_constraint(real mu, real rho, real nu, int N){
        vector[N+1] lp ;
        real h ;
        h = 2 * pi() / N ;
        lp[1] = aewrappedcauchy_lpdf(-pi() | mu, rho, nu) ;
        for (n in 1:(N/2)){
            lp[2*n] = log(4) + aewrappedcauchy_lpdf(-pi() + h*(2*n-1) | mu, rho, nu) ;
        }
        for (n in 1:(N/2-1)){
            lp[2*n+1] = log(2) + aewrappedcauchy_lpdf(-pi() + h*2*n | mu, rho, nu) ;
        }
        lp[N+1] = aewrappedcauchy_lpdf(pi() | mu, rho, nu) ;
        return (log(h/3) + log_sum_exp(lp)) ;

    } 

    real aewrappedcauchy_mixture_lpdf(real R, int K, vector a, vector mu, vector rho, vector nu) {
        vector[K] lp ;
        for (k in 1:K){
            real logncon ;
            logncon = aevon_mises_normalize_constraint(mu[k], rho[k], nu[k], 20) ;
            lp[k] = log(a[k]) + aewrappedcauchy_lpdf(R | mu[k], rho[k], nu[k]) - logncon ;
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
        rho[s] ~ student_t(2.5, 0, 0.17) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * aewrappedcauchy_mixture_lpdf(RADIAN[i] | K, alpha, ori, rho[SUBJECT[i]], nu[SUBJECT[i]]) ;
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
        log_lik[i] = DEPTH[i] * aewrappedcauchy_mixture_lpdf(RADIAN[i] | K, alpha, ori, rho[SUBJECT[i]], nu[SUBJECT[i]]) ;
    }
}
