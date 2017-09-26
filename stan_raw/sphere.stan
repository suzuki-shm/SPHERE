data {
    int I ;
    int D[I] ;
}

parameters {
    real y0 ;
    real<lower=0> H ;
    real<lower=-1, upper=1> O[2] ;
    vector<lower=-pi()/2, upper=pi()/2>[I-1] y_raw ;
    real<lower=0> sigma ;
}

transformed parameters{
    vector[I] y ;
    vector[I] lambda ;
    vector[I] trend ;

    // variance
    y[1] = y0 ;
    for(i in 2:I){
        y[i] = y[i-1] + sigma*tan(y_raw[i-1]) ;
    }
    // trend from replication rate
    for(i in 1:I){
        trend[i] = H * (cos((i-atan2(O[1], O[2])/(2*pi())*I)*2*pi()/I) + 1) ;
    }
    lambda = exp(y + trend) ;

}

model {
    D ~ poisson(lambda) ;
}

generated quantities {
    real PTR ;
    PTR = exp(H*2) ;
}
