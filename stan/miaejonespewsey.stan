functions {
    real trans_sin2(real theta, real mu, real nu){
        real theta_mu ;
        theta_mu =  theta - nu * sin(theta - mu) * sin(theta - mu) ;
        return theta_mu ;
    }

    real jonespewsey_lpdf(real theta, real mu, real kappa, real psi){
        if (fabs(psi) < 1e-10){
            return kappa * cos(theta - mu) ;
        }else{
            return 1 / psi * log(cosh(kappa * psi) + sinh(kappa * psi) * cos(theta - mu)) ;
        }
    }

    real miaejonespewsey_lpdf(real theta, real mu, real kappa, real psi, real nu){
        return jonespewsey_lpdf(trans_sin2(theta, mu, nu)| mu, kappa, psi) ;
    }

    real miaejonespewsey_normalize_constraint(real mu, real kappa, real psi, real nu, int N){
        // Numerical integration by composite Simpson's rule
        vector[N+1] lp ;
        real h ;

        h = 2 * pi() / N ;
        lp[1] = miaejonespewsey_lpdf(-pi() | mu, kappa, psi, nu) ;
        for (n in 1:(N/2)){
            lp[2*n] = log(4) + miaejonespewsey_lpdf(-pi() + h*(2*n-1) | mu, kappa, psi, nu) ;
        }
        for (n in 1:(N/2-1)){
            lp[2*n+1] = log(2) + miaejonespewsey_lpdf(-pi() + h*2*n | mu, kappa, psi, nu) ;
        }
        lp[N+1] = miaejonespewsey_lpdf(pi() | mu, kappa, psi, nu) ;
        return log(h/3) + log_sum_exp(lp) ;
    }

    real miaejonespewsey_mixture_lpdf(real R, int K, vector a, vector mu, vector kappa, vector psi, vector nu){
        vector[K] lp ;
        real logncon ;

        for (k in 1:K){
            logncon = miaejonespewsey_normalize_constraint(mu[k], kappa[k], psi[k], nu[k], 20) ;
            lp[k] = log(a[k]) + miaejonespewsey_lpdf(R | mu[k], kappa[k], psi[k], nu[k]) - logncon ;
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
    vector<lower=-1.0, upper=1.0>[K] psi[S] ;
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
        kappa[s] ~ student_t(2.5, 0, 0.2) ;
    }
    for(i in 1:I){
        target += DEPTH[i] * miaejonespewsey_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], psi[SUBJECT[i]], nu[SUBJECT[i]]) ;
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
        log_lik[i] = DEPTH[i] * miaejonespewsey_mixture_lpdf(RADIAN[i] | K, alpha, ori, kappa[SUBJECT[i]], psi[SUBJECT[i]], nu[SUBJECT[i]]) ;
    }
}
