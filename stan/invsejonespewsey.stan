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
                    t1 = t ;
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

    real jonespewsey_lpdf(real theta, real mu, real kappa, real psi){
        if (fabs(psi) < 1e-10){
            return kappa * cos(theta - mu) ;
        }else{
            return 1.0 / psi * log(cosh(kappa * psi) + sinh(kappa * psi)*cos(theta - mu)) ;
        }
    }

    real invsejonespewsey_lpdf(real theta, real mu, real kappa, real psi, real lambda){
        return jonespewsey_lpdf(trans_t_lambda(theta, lambda, mu) | mu, kappa, psi) ;
    }

    real g_jp(real theta, real kappa, real psi, real lambda){
        return log(1.0 - (1.0 + lambda) * cos(theta) / 2.0) + jonespewsey_lpdf(theta - (1.0 - lambda) * sin(theta) / 2.0| 0, kappa, psi) ;
    }

    real g_vm(real theta, real kappa, real lambda){
        return log(1.0 - (1.0 + lambda) * cos(theta) / 2.0) + von_mises_lpdf(theta - (1.0 - lambda) * sin(theta) / 2.0| 0, kappa) ;
    }

    real invsejonespewsey_normalize_constraint(real mu, real kappa, real psi, real lambda, int N){
        // Numerical integration by composite Simpson's rule
        vector[N+1] lp ;
        real h ;
        real logncon ;

        h = 2.0 * pi() / N ;
        if (fabs(psi) < 1e-10){
            lp[1] = g_vm(-pi(), kappa, lambda) ;
            for (n in 1:N/2){
                lp[2*n] = log(4) + g_vm(-pi() + h*(2*n-1), kappa, lambda) ;
            }
            for (n in 1:N/2-1){
                lp[2*n+1] = log(2) + g_vm(-pi() + h*2*n, kappa, lambda) ;
            }
            lp[N+1] = g_vm(pi(), kappa, lambda) ;
        }else{
            lp[1] = g_jp(-pi(), kappa, psi, lambda) ;
            for (n in 1:N/2){
                lp[2*n] = log(4) + g_jp(-pi() + h*(2*n-1), kappa, psi, lambda) ;
            }
            for (n in 1:N/2-1){
                lp[2*n+1] = log(2) + g_jp(-pi() + h*2*n, kappa, psi, lambda) ;
            }
            lp[N+1] = g_jp(pi(), kappa, psi, lambda) ;
        }
        logncon =  log(h/3.0) + log_sum_exp(lp) ;
        return logncon ;
    }

    real invsejonespewsey_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, vector psi, vector lambda){
        vector[K] lp ;
        real logncon ;

        for (k in 1:K){
            logncon = invsejonespewsey_normalize_constraint(mu[k], kappa[k], psi[k], lambda[k], 20) ;
            lp[k] = log(a[k]) + invsejonespewsey_lpdf(R | mu[k], kappa[k], psi[k], lambda[k]) - logncon ;
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

    RADIAN = -pi() + (2.0 * pi() / L) * (to_vector(LOCATION) - 1.0) ;
    A = rep_vector(50.0/K, K) ;
}

parameters {
    simplex[K] alpha ;
    unit_vector[2] O[K] ;
    vector<lower=0.0>[K] kappa[S] ;
    vector<lower=-1.0, upper=1.0>[K] psi[S] ;
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
        psi[s] ~ normal(0, 1.0) ;
        lambda[s] ~ normal(0, 1.0) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * invsejonespewsey_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], psi[SUBJECT[i]], lambda[SUBJECT[i]]) ;
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
        log_lik[i] = DEPTH[i] * invsejonespewsey_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], psi[SUBJECT[i]], lambda[SUBJECT[i]]) ;
    }
}
