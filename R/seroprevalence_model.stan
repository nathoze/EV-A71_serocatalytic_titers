data {
    int <lower=0> NFOI; //the number of FOI estimated (30: 1983 - 2012)   
    int<lower = 0> NAges; //Number of age classes (12: 1 - 12)
    int<lower = 0> Nyears; //Number of years observed (18: 1995 - 2012)

}

parameters {
  real<lower=0> FOI;
  real<lower=0> risk_H;
}


transformed parameters {
    real<lower =0, upper=1> pinf;
//   real<lower =0, upper=1> P;
    real<lower =0, upper=1> P[Nyears,NAges];

    pinf = 1-exp(-FOI[j]);

    P = 0
    for(j = 0; j < Nyears; j++){
       for(i = 0; i< NAges; i++){
         pinf = 1-exp(-FOI[j]);
          P[j,i] = pinf;
      }    
    }
 
}
model {
    
    for (j in 1:NFOI) {    
      FOI[j] ~ exponential(2);
    }
  omega ~ exponential(2);
  
    target += xxx

}

