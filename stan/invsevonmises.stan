functions {
    real inv_trans_batschelet(real theta, real mu, real lambda){
        real t ;
        real ft ;
        real err ;
        int count ;
        count = 0 ;
        // Small lambda works with Newton's method
        if (fabs(lambda) <= 0.8){
            t = theta ;
            ft = t - (1+lambda) * sin(t-mu) / 2 - theta ;
            err = fabs(ft) ;
            while(err > 1e-8){
                t = t - ( ft / (1 - (1+lambda) * cos(t-mu) / 2)) ;
                ft = t - (1+lambda) * sin(t-mu) / 2 - theta ;
                err = fabs(ft) ;
                count += 1 ;
                if (count >= 30){
                    break ;
                }
            }
        // Large lambda only works with bisection method
        }else{
            real t1 ;
            real t2 ;
            t1 = -2.0*pi() ;
            t2 = 2.0*pi() ;
            t = (t1 + t2) / 2 ;
            ft = t - (1+lambda) * sin(t-mu) / 2 - theta  ;
            err = fabs(ft) ;
            while(err > 1e-8){
                if (ft < 0){
                    t1 = t ;
                }else{
                    t2 = t ;
                }
                t = (t1 + t2) / 2 ;
                ft = t - (1+lambda) * sin(t-mu) / 2 - theta  ;
                err = fabs(ft) ;
                count += 1 ;
                if (count >= 50){
                    break ;
                }
            }
        }
        return t ;
    }

    real trans_t_lambda(real theta, real lambda, real mu){
        real t_lambda ;
        t_lambda = (1 - lambda) / (1 + lambda) * theta + 2 * lambda / (1 + lambda) * inv_trans_batschelet(theta, mu, lambda) ;
        return t_lambda ;
    }

    real invsevon_mises_lpdf(real theta, real mu, real kappa, real lambda){
        return von_mises_lpdf(trans_t_lambda(theta, lambda, mu) | mu, kappa) ;
    }

    real g(real theta, real kappa, real lambda){
        return log(1 - (1 + lambda) * cos(theta) / 2) + von_mises_lpdf(theta - (1 - lambda) * sin(theta) / 2| 0, kappa) ;
    }

    real invsevon_mises_normalize_constraint(real kappa, real lambda, int N){
        // Numerical integration by composite Simpson's rule
        vector[N+1] lp ;
        real h ;

        h = 2 * pi() / N ;
        lp[1] = g(-pi(), kappa, lambda) ;
        for (n in 1:N/2){
            lp[2*n] = log(4) + g(-pi() + h*(2*n-1), kappa, lambda) ;
        }
        for (n in 1:N/2-1){
            lp[2*n+1] = log(2) + g(-pi() + h*2*n, kappa, lambda) ;
        }
        lp[N+1] = g(pi(), kappa, lambda) ;
        return (log(h/3) + log_sum_exp(lp)) ;

    }

    real invsevon_mises_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, vector lambda){
        vector[K] lp;
        real logncon ;
        for (k in 1:K){
            logncon = invsevon_mises_normalize_constraint(kappa[k], lambda[k], 20) ;
            lp[k] = log(a[k]) + invsevon_mises_lpdf(R | mu[k], kappa[k], lambda[k]) - logncon ;
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
    vector<lower=0.0>[K] kappa[S] ;
    // peakedness parameter
    vector<lower=-1.0, upper=1.0>[K] lambda[S] ;
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
        kappa[s] ~ student_t(2.5, 0, 0.2025) ;
        lambda[s] ~ normal(0, 1.0) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * invsevon_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], lambda[SUBJECT[i]]) ;
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
        log_lik[i] = DEPTH[i] * invsevon_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
}
