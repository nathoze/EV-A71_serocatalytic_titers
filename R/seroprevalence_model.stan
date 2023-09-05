data {
  int <lower=0> NFOI; //the number of FOI estimated (30: 1983 - 2012)
  int<lower = 0> NAges; //Number of age classes (12: 1 - 12)
  int<lower = 0> Nyears; //Number of years observed (18: 1995 - 2012)
  int<lower =0> nsamples[Nyears,NAges];
  int<lower =0> seropositive[Nyears,NAges];
}

parameters {
  real<lower=0> FOI[NFOI];
  real<lower=0> omega;
}


transformed parameters {
  real<lower =0, upper=1> pinf;
  real<lower =0, upper=1> P[Nyears,NAges];
  //  array[Nyears,NAges] real P;

  for(birthyear in 1:NFOI){
    pinf = 0;
    for(age in 0:NAges){
    //  sampling year  = i+birthyear-1
      if(age+birthyear <=NFOI){
        pinf=pinf*(exp(-FOI[age+birthyear]-omega)) + FOI[age+birthyear]/(FOI[age+birthyear]+omega)*(1-exp(-omega-FOI[age+birthyear]));
        if(age >0 && (age+birthyear-12>= 1)){
          P[age+birthyear-12,age] = pinf;
        }
      }
    }
  }

}

model {

  for (j in 1:NFOI) {
    FOI[j] ~ exponential(2);
  }
  omega ~ exponential(2);
  for(j in 1:Nyears){
    for(i in 1:NAges){
      target += binomial_lpmf(seropositive[j,i]|nsamples[j,i], P[j,i]) ;
    }
  }

}

