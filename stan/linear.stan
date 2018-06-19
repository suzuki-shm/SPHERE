data {
    int I ;
    int S ;
    int L ;
    int<lower=1, upper=L> LOCATION[I] ;
    int<lower=1, upper=S> SUBJECT[I] ;
    int<lower=0> DEPTH[I] ;
}

parameters {
    unit_vector[2] O ;
    vector<lower=0.0>[S] H ;
    real<lower=0.0> sigma_H ;
    real flex ;
}

transformed parameters{
    real<lower=-pi(), upper=pi()> ori ;
    vector[L] trend[S] ;
    vector<lower=0.0>[L] lambda[S] ;

    // convert unit vector
    ori = atan2(O[1], O[2]) ;
    for(s in 1:S){
        // trend from replication rate
        for(l in 1:L){
            trend[s, l] = 2.0 * H[s] / L * fabs(fabs(l - ori / 2.0 / pi() * L) - L / 2.0) ;
        }
        lambda[s] = exp(trend[s] + flex);
    }
}

model {
    H ~ normal(0, sigma_H) ;
    for(i in 1:I){
        DEPTH[i] ~ poisson(lambda[SUBJECT[i], LOCATION[i]]) ;
    }
}

generated quantities {
    vector<lower=1.0>[K] PTR[S] ;
    vector[I] log_lik ;

    PTR = exp(H) ;
    for(i in 1:I){
        log_lik[i] = poisson_lpmf(DEPTH[i] | lambda[SUBJECT[i], LOCATION[i]]) ;
    }
}
