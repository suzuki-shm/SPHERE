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
            ft = t - (1.0+lambda) * sin(t-mu) / 2.0 - theta ;
            err = fabs(ft) ;
            while(err > machine_precision()){
                t = t - ( ft / (1.0 - (1.0+lambda) * cos(t-mu) / 2.0)) ;
                ft = t - (1.0+lambda) * sin(t-mu) / 2.0 - theta ;
                err = fabs(ft) ;
                count += 1 ;
                if (count >= 30){
                    break ;
                }
            }
        // Large lambda only works with illinois method
        }else{
            real t1 ;
            real t2 ;
            real ft1 ;
            real ft2 ;
            t1 = -2.0*pi() ;
            t2 = 2.0*pi() ;
            ft1 = t1 - (1.0+lambda) * sin(t1-mu) / 2.0 - theta  ;
            ft2 = t2 - (1.0+lambda) * sin(t2-mu) / 2.0 - theta  ;
            t = (t1 * ft2 - t2 * ft1) / (ft2 - ft1) ;
            ft = t - (1.0+lambda) * sin(t-mu) / 2.0 - theta  ;
            err = fabs(ft) ;
            // By loss of significant digits, this method cannot give lower error than machine precision
            while(err > 1e-13){
                if (ft < 0){
                    t1 = t ;
                    ft1 = t1 - (1.0+lambda) * sin(t1-mu) / 2.0 - theta  ;
                    ft2 /= 2.0 ;
                }else{
                    t2 = t ;
                    ft1 /= 2.0 ;
                    ft2 = t2 - (1.0+lambda) * sin(t2-mu) / 2.0 - theta  ;
                }
                t = (t1 * ft2 - t2 * ft1) / (ft2 - ft1) ;
                ft = t - (1.0+lambda) * sin(t-mu) / 2.0 - theta  ;
                err = fabs(ft) ;
                count += 1 ;
                if (count >= 100){
                    break ;
                }
            }
        }
        return t ;
    }

    real trans_t_lambda(real theta, real lambda, real mu){
        real t_lambda ;
        t_lambda = (1.0 - lambda) / (1.0 + lambda) * theta + 2.0 * lambda / (1.0 + lambda) * inv_trans_batschelet(theta, mu, lambda) ;
        return t_lambda ;
    }

    real invsevon_mises_lpdf(real theta, real mu, real kappa, real lambda){
        return von_mises_lpdf(trans_t_lambda(theta, lambda, mu) | mu, kappa) ;
    }

    real g(real theta, real kappa, real lambda){
        return log(1 - (1 + lambda) * cos(theta) / 2.0) + von_mises_lpdf(theta - (1 - lambda) * sin(theta) / 2.0| 0, kappa) ;
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
        return (log(h/3.0) + log_sum_exp(lp)) ;

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
    // Unconstrained concentration parameter
    vector[K] kappa_uncon[S] ;
    // peakedness parameter
    vector<lower=-1.0, upper=1.0>[K] lambda[S] ;
}

transformed parameters{
    vector[K] ori ;
    vector<lower=0.0, upper=4.0>[K] kappa[S] ;

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
        lambda[s] ~ normal(0, 1.0) ;
        // Jacobian adjustment for parameter transformation (see 'Lower and Upper Bounded Scalar' in Stan manual)
        target += log(log(4) ./ (2 * alpha)) + log_inv_logit(kappa_uncon[s]) + log1m_inv_logit(kappa_uncon[s]) ;
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
        PTR[s] = exp(2.0 * kappa[s]) ;
        wPTR[s] = exp(2.0 * alpha .* kappa[s]) ;
        mwPTR[s] = sum(wPTR[s]) ;
    }
    for(i in 1:I){
        log_lik[i] = DEPTH[i] * invsevon_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
}
