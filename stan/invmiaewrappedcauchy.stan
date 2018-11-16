functions{
    real inv_trans_sin2(real theta, real mu, real nu){
        real t ;
        real ft ;
        real err ;
        int count ;
        count = 0 ;
        // Small nu works with Newton's method
        if (fabs(nu) <= 0.8){
            t = theta ;
            ft = t + nu*pow(sin(t-mu),2) - theta ;
            err = fabs(ft) ;
            while(err > 1e-8){
                t = t - (ft / (1 + 2 * nu * sin(t-mu) * cos(t-mu))) ;
                ft = t + nu*pow(sin(t-mu),2) - theta ;
                err = fabs(ft) ;
                count += 1 ;
                if (count >= 30){
                    break ;
                }
            }
        // Large nu only works with bisection method
        }else{
            real t1 ;
            real t2 ;
            t1 = -2.0*pi() ;
            t2 = 2.0*pi() ;
            t = (-pi() + pi()) / 2 ;
            ft = t + nu*pow(sin(t-mu),2) - theta ;
            err = fabs(ft) ;
            while(err > 1e-8){
                if (ft < 0){
                    t1 = t ;
                }else{
                    t2 = t ;
                }
                t = (t1 + t2) / 2 ;
                ft = t + nu*pow(sin(t-mu),2) - theta ;
                err = fabs(ft) ;
                count += 1 ;
                if (count >= 50){
                    break ;
                }
            }
        }
        return t ;
    }

    real wrappedcauchy_lpdf(real theta, real mu, real rho){
        return log(1 - pow(rho, 2)) - log(2) - log(pi()) - log(1 + pow(rho, 2) - 2 * rho * cos(theta - mu)) ;
    }

    real invmiaewrappedcauchy_mixture_lpdf(real R, int K, vector a, vector mu, vector rho, vector nu) {
        vector[K] lp;
        for (k in 1:K){
            lp[k] = log(a[k]) + wrappedcauchy_lpdf(inv_trans_sin2(R, mu[k], nu[k]) | mu[k], rho[k]) ;
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
        rho[s] ~ student_t(2.5, 0, 0.10125) ;
        nu[s] ~ normal(0, 1.0) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * invmiaewrappedcauchy_mixture_lpdf(RADIAN[i] | K, alpha, ori, rho[SUBJECT[i]], nu[SUBJECT[i]]) ;
    }
}

generated quantities {
    vector<lower=0.0>[K] kappa[S] ;
    vector<lower=1.0>[K] PTR[S] ;
    vector<lower=1.0>[K] wPTR[S] ;
    vector<lower=1.0>[S] mwPTR ;
    vector[I] log_lik ;

    for(s in 1:S){
        // See (Jones&Pewsey, 2005) about this transformation
        kappa[s] = 2 * atanh(rho[s]) ;
        // Fold change of max p.d.f. to min p.d.f.
        PTR[s] = exp(2 * kappa[s]) ;
        wPTR[s] = exp(2 * alpha .* kappa[s]) ;
        mwPTR[s] = sum(wPTR[s]) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * invmiaewrappedcauchy_mixture_lpdf(RADIAN[i] | K, alpha, ori, rho[SUBJECT[i]], nu[SUBJECT[i]]) ;
    }
}
