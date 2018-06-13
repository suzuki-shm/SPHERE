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
    vector<lower=0>[S] H ;
    real flex0[S] ;
    vector<lower=-pi()/2, upper=pi()/2>[L-1] flex_raw[S] ;
    real<lower=0> sigma_flex[S] ;
    real<lower=0> sigma_sigma_flex ;
    real<lower=0> sigma_H ;
}

transformed parameters{
    real<lower=-pi(), upper=pi()> ori ;
    vector[L] flex[S] ;
    vector[L] trend[S] ;
    vector<lower=0>[L] lambda[S] ;
    // convert unit vector
    ori = atan2(O[1], O[2]) ;
    for(s in 1:S){
        // flex
        flex[s, 1] = flex0[s] ;
        for(l in 2:L){
            flex[s, l] = flex[s, l-1] + sigma_flex[s] * tan(flex_raw[s, l-1]) ;
        }
        // trend from replication rate
        for(l in 1:L){
            trend[s, l] = H[s] / 2.0 * (cos(l * 2.0 * pi() / L - ori) + 1.0) ;
        }
        lambda[s] = exp(flex[s] + trend[s]) ;
    }
}

model {
    H ~ normal(0, sigma_H) ;
    for(s in 1:S){
        sigma_flex[s] ~ normal(0, sigma_sigma_flex) ;
    }
    for(i in 1:I){
        DEPTH[i] ~ poisson(lambda[SUBJECT[i], LOCATION[i]]) ;
    }
}

generated quantities {
    vector<lower=1.0>[S] PTR ;
    vector[I] log_lik ;
    PTR = exp(H) ;
    for(i in 1:I){
        log_lik[i] = poisson_lpmf(DEPTH[i] | lambda[SUBJECT[i], LOCATION[i]]) ;
    }
}
