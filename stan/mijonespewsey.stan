functions{
    real jonespewsey_lpdf(real theta, real mu, real kappa, real psi){
        if (fabs(psi) < 1e-10){
            return kappa * cos(theta - mu) ;
        }else{
            return 1 / psi * log(cosh(kappa * psi) + sinh(kappa * psi)*cos(theta - mu)) ;
        }
    }

    real jonespewsey_normalize_constraint(real mu, real kappa, real psi, int N){
        // Numerical integration by composite Simpson's rule
        vector[N+1] lp ;
        real h ;
        real logncon ;

        if (fabs(psi) < 1e-10){
            logncon = log(2 * pi() * modified_bessel_first_kind(0, kappa)) ;
        }else{
            h = 2 * pi() / N ;
            lp[1] = jonespewsey_lpdf(-pi() | mu, kappa, psi) ;
            for (n in 1:(N/2)){
                lp[2*n] = log(4) + jonespewsey_lpdf(-pi() + h*(2*n-1) | mu, kappa, psi) ;
            }
            for (n in 1:(N/2-1)){
                lp[2*n+1] = log(2) + jonespewsey_lpdf(-pi() + h*2*n | mu, kappa, psi) ;
            }
            lp[N+1] = jonespewsey_lpdf(pi() | mu, kappa, psi) ;
            logncon =  log(h/3) + log_sum_exp(lp) ;
        }
        return logncon ;
    }

    // system of inverse Batschelet transformation
    vector invAPF(vector y_init, vector t, real[] x_r, int[] x_i){
        // Get length of array
        int I = x_i[1] ;
        vector[I] y ; // value to be zero
        vector[I] theta = t[:I] ; // known angle
        vector[I] nu = t[I+1:] ;

        // construct equation to be zero
        y = -theta + y_init + nu .* sin(y_init) .* sin(y_init);
        return y ;
    }
}

data {
    int I ;
    int S ;
    int L ;
    int<lower=1, upper=S> SUBJECT[I] ;
    int<lower=1, upper=L> LOCATION[I] ;
    int<lower=0> DEPTH[I] ;
    int<lower=1> K ;
}

transformed data {
    real x_r[0] ;
    int x_i[1] ;
    vector<lower=-pi(), upper=pi()>[I] RADIAN ;
    vector<lower=-pi(), upper=pi()>[I] THETA_INIT ;
    vector<lower=0.0>[K] A;

    // set integer parameter to pass solver
    x_i[1] = I ;
    // set initial value of inverse transformed value
    for (i in 1:I){
        RADIAN[i] = -pi() + (2.0 * pi() / L) * (LOCATION[i] - 1) ;
        THETA_INIT[i] = 0 ;
    }
    for(k in 1:K){
        A[k] = 50 / K ;
    }
}

parameters {
    simplex[K] alpha ;
    unit_vector[2] O[K] ;
    vector<lower=0.0>[K] kappa[S] ;
    vector<lower=-1.0, upper=1.0>[K] psi[S] ;
    vector<lower=-1.0, upper=1.0>[K] nu[S] ;
}

transformed parameters{
    vector[K] ori ;
    vector[I] THETA[K] ;
    vector[2*I] t ;

    for (k in 1:K){
        ori[k] = atan2(O[k][1], O[k][2]) ;
        // construct vector for algebra solver
        for (i in 1:I){
            // transform angle by location parameter
            t[i] = RADIAN[i] - ori[k] ;
            if (t[i] > pi()){
                t[i] = t[i] - 2*pi() ;
            }else if (t[i] <= -pi()){
                t[i] = t[i] + 2*pi() ;
            }
            t[I+i] = nu[SUBJECT[i]][k] ;
        }
        THETA[k] = algebra_solver(invAPF, THETA_INIT, t, x_r, x_i) ;
    }
}

model {
    alpha ~ dirichlet(A) ;
    for(s in 1:S){
        kappa[s] ~ student_t(2.5, 0, 0.2) ;
    }
    for (i in 1:I){
        vector[K] lp;
        real logncon ;
        for (k in 1:K){
            logncon = jonespewsey_normalize_constraint(ori[k], kappa[SUBJECT[i]][k], psi[SUBJECT[i]][k], 20) ;
            lp[k] = log(alpha[k]) + jonespewsey_lpdf(THETA[k][i] | 0, kappa[SUBJECT[i]][k], psi[SUBJECT[i]][k]) -logncon ;
        }
        target += DEPTH[i] * log_sum_exp(lp) ;
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
        vector[K] lp;
        real logncon ;
        for (k in 1:K){
            logncon = jonespewsey_normalize_constraint(ori[k], kappa[SUBJECT[i]][k], psi[SUBJECT[i]][k], 20) ;
            lp[k] = log(alpha[k]) + jonespewsey_lpdf(THETA[k][i] | 0, kappa[SUBJECT[i]][k], psi[SUBJECT[i]][k]) -logncon ;
        }
        log_lik[i] = DEPTH[i] * log_sum_exp(lp) ;
    }
}

