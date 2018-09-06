functions{
    // system of inverse Abe-Pewsey-Fujisawa transformation
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
    for (i in 1:I){
        RADIAN[i] = -pi() + (2.0 * pi() / L) * (LOCATION[i] - 1) ;
        // set initial value of inverse transformed value
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
        for (k in 1:K){
            lp[k] = log(alpha[k]) + von_mises_lpdf(THETA[k][i] | 0, kappa[SUBJECT[i]][k]) ;
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
        for (k in 1:K){
            lp[k] = log(alpha[k]) + von_mises_lpdf(THETA[k][i] | 0, kappa[SUBJECT[i]][k]) ;
        }
        log_lik[i] = DEPTH[i] * log_sum_exp(lp) ;
    }
}

