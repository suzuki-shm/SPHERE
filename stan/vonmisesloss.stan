functions {
    real von_mises_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa) {
        vector[K] lp;
        for (k in 1:K){
            lp[k] = log(a[k]) + von_mises_lpdf(R | mu[k], kappa[k]) ;
        }
        return log_sum_exp(lp) ;
    }

    real uniform_loglik_Simpson(real lower, real upper){
        return log((upper - lower) / 2 / pi()) ;
    }

    real von_mises_loglik_Simpson(real lower, real upper, int K, vector a, vector mu, vector kappa){
        int M = 20;
        vector[M+1] lp;
        real h;
        h = (upper - lower) / M ;
        lp[1] = von_mises_mixture_lpdf(lower | K, a, mu, kappa) ;
        for (m in 1:M/2){
            lp[2*m] = log(4) + von_mises_mixture_lpdf(lower + h*(2*m-1) | K, a, mu, kappa) ;
        }
        for (m in 1:M/2-1){
            lp[2*m+1] = log(2) + von_mises_mixture_lpdf(lower + h*2*m | K, a, mu, kappa) ;
        }
        lp[M+1] = von_mises_mixture_lpdf(upper | K, a, mu, kappa) ;
        return (log(h/3) + log_sum_exp(lp)) ;
    }
}

data {
    int I ; // Input length
    int L ; // Complete genome length
    int S ; // Subject number
    int C ; // Contig number
    int<lower=1, upper=L> L_contig[C] ; // Contig length
    int<lower=1, upper=S> SUBJECT[I] ;
    int<lower=1> CONTIG[I] ;
    int<lower=1> LOCATION[I] ;
    int<lower=0> DEPTH[I] ;
    int<lower=1> K ; // number of mixed distribution
}

transformed data {
    vector<lower=-pi(), upper=pi()>[I] RADIAN ;
    vector<lower=0.0>[K] A; //hyperparameter for dirichlet distribution
    int DEPTH_SUM[S, C] ;
    // Integrate of uniform distribution in the range
    real S_uniform[C] ;

    for (i in 1:I){
        RADIAN[i] = -pi() + (2.0 * pi() / L) * (LOCATION[i] - 1) ;
    }

    A = rep_vector(50.0/K, K) ;
    for (s in 1:S){
        for (c in 1:C){
            DEPTH_SUM[s, c] = 0 ;
        }
    }
    for (i in 1:I){
        DEPTH_SUM[SUBJECT[i], CONTIG[i]] += DEPTH[i] ;
    }
    for (c in 1:C){
        S_uniform[c] = uniform_loglik_Simpson(-pi(), -pi() + (2.0 * pi() / L) * (L_contig[c] - 1)) ;
    }
}

parameters {
    simplex[K] alpha ;
    unit_vector[2] O[C, K] ;
    // Unconstrained concentration parameter
    vector[K] kappa_uncon[S] ;
    real<lower=0> sigma[S] ;
}

transformed parameters{
    vector[K] ori[C] ;
    vector<lower=0>[K] kappa[S] ;
    real S_von_mises[S, C] ;

    for (c in 1:C){
        for (k in 1:K){
            // convert unit vector
            ori[c][k] = atan2(O[c, k][1], O[c, k][2]) ;
        }
    }

    for (s in 1:S){
        // Add upper bound to kappa using alpha (see 'Lower and Upper Bounded Scalar' in Stan manual)
        kappa[s] = log(4) ./ (2 * alpha) .* inv_logit(kappa_uncon[s]) ;
        // Integrate of circular distribution in the range
        for (c in 1:C){
            S_von_mises[s, c] = von_mises_loglik_Simpson(-pi(), -pi() + (2.0 * pi() / L) * (L_contig[c] - 1), K, alpha, ori[c], kappa[s]) ;
        }
    }
}

model {    
    alpha ~ dirichlet(A) ;
    for (s in 1:S){
        alpha .* kappa[s] ~ student_t(2.5, 0, 0.2025) ;
        // Jacobian adjustment for parameter transformation (see 'Lower and Upper Bounded Scalar' in Stan manual)
        target += log(log(4) ./ (2 * alpha)) + log_inv_logit(kappa_uncon[s]) + log1m_inv_logit(kappa_uncon[s]) ;
        for (c in 1:C){
            // Prior distribution for sigma
            sigma[s] ~ student_t(2.5, 0, DEPTH_SUM[s, c] * (L_contig[c]-1) / pow(L_contig[c], 2)) ;
        }
    }
    for(i in 1:I){
        DEPTH[i] ~ normal(DEPTH_SUM[SUBJECT[i], CONTIG[i]] * exp(von_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori[CONTIG[i]], kappa[SUBJECT[i]]) + log(2 * pi()) + S_uniform[CONTIG[i]] - log(L_contig[CONTIG[i]]) - S_von_mises[SUBJECT[i], CONTIG[i]]), sigma[SUBJECT[i]]) ;
    }
}

generated quantities {
    vector<lower=1.0>[K] PTR[S] ;
    vector<lower=1.0>[K] wPTR[S] ;
    vector<lower=1.0>[S] mwPTR ;
    vector<lower=0.0, upper=1.0>[K] MRL[S] ;
    vector<lower=0.0, upper=1.0>[K] CV[S] ;
    vector<lower=0.0>[K] CSD[S] ;
    vector[I] log_lik ;

    for(s in 1:S){
        // Fold change of max p.d.f. to min p.d.f.
        PTR[s] = exp(2 * kappa[s]) ;
        wPTR[s] = exp(2 * alpha .* kappa[s]) ;
        mwPTR[s] = sum(wPTR[s]) ;
        // Mean resultant length
        for (k in 1:K){
            MRL[s][k] = modified_bessel_first_kind(1, kappa[s][k]) / modified_bessel_first_kind(0, kappa[s][k]) ;
        }
        // Circular variance
        CV[s] = 1 - MRL[s] ;
        // Circular standard variation
        CSD[s] = sqrt(-2 * log(MRL[s])) ;
    }
    for(i in 1:I){
        log_lik[i] = normal_lpdf(DEPTH[i] | DEPTH_SUM[SUBJECT[i], CONTIG[i]] * exp(von_mises_mixture_lpdf(RADIAN[i] | K, alpha, ori[CONTIG[i]], kappa[SUBJECT[i]]) + log(2 * pi()) + S_uniform[CONTIG[i]] - log(L_contig[CONTIG[i]]) - S_von_mises[SUBJECT[i], CONTIG[i]]), sigma[SUBJECT[i]]) ;
    }
}
