functions {
    real jonespewsey_lpdf(real theta, real mu, real kappa, real psi){
        if (fabs(psi) < 1e-10){
            return kappa * cos(theta - mu) ;
        }else{
            return 1 / psi * log(cosh(kappa * psi) + sinh(kappa * psi)*cos(theta - mu)) ;
        }
    }

    real djonespewsey_normalize_constraint(real mu, real kappa, real psi, int N){
        vector[N] lp;

        for (n in 1:N){
            real theta ;
            theta = -pi() + (2.0 * pi() / N) * (n - 1) ;
            lp[n] = jonespewsey_lpdf(theta | mu, kappa, psi) ;
        }
        return log_sum_exp(lp) ;
    }

    real djonespewsey_lpdf(real theta, real mu, real kappa, real psi, int N){
        real logncon ;
        logncon = djonespewsey_normalize_constraint(mu, kappa, psi, N) ;
        return jonespewsey_lpdf(theta | mu, kappa, psi) - logncon ;
    }

    real djonespewsey_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, vector psi, int N){
        vector[K] lp ;
        for (k in 1:K){
            lp[k] = log(a[k]) + djonespewsey_lpdf(R | mu[k], kappa[k], psi[k], N) ;
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
    vector<lower=0.0>[K] kappa[S] ;
    vector<lower=-1.0, upper=1.0>[K] psi[S] ;
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
        target += DEPTH[i] * djonespewsey_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], psi[SUBJECT[i]], L) ;
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
        log_lik[i] = DEPTH[i] * djonespewsey_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], psi[SUBJECT[i]], L) ;
    }
}
