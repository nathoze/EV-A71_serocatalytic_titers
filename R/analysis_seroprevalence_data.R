rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(purrr)
library(Hmisc)

source("R/Model_EVA71.R")
source("R/utils.R")
source("R/FOI_models.R")
source("R/plots.R")

################################################################################
# Data : EV71 in Malaysia 1995 -> 2012
################################################################################

data.EV71.Malaysia  <- readRDS('data/titers_by_class.rds')
data.EV71.Malaysia <- data.EV71.Malaysia %>%
  mutate(age = sampling.year-birth.year) %>%
  mutate(n =  replace_na(n,0) ) %>%
  filter(age <= 12)

################################################################################
# Parameters
################################################################################

min.year.sampling <- min(data.EV71.Malaysia$sampling.year)
max.year.sampling <- max(data.EV71.Malaysia$sampling.year)
sampling.years = seq(min.year.sampling, max.year.sampling)
N.sampling.years = length(sampling.years)
age.min = 1 # minimal age in all the surveys at the time of survey
age.max = 12 # maximal age in all surveys at the time of survey
N.FOI = N.sampling.years + age.max-age.min # The number of FOI measured. We consider here that it is possible to be infected at age 0
N.birth.years   = N.FOI # The number of birth years considered (-1 if we remove the latest year)
birth.years = seq(min.year.sampling-age.max, min.year.sampling-age.max+N.FOI-1)

#start.sampling.year = age.max
N.titer.sets = age.max*N.sampling.years
titer.observable.max = max(data.EV71.Malaysia$titer.class)

N.titers = 10 # the number of possible values of the titers
Titers.0 <- c(1,rep(0,N.titers-1)) # probability distribution of the titers at birth (everybody is seronegative at birth)
data.parameters = list(min.year.sampling = min.year.sampling,
                       max.year.sampling = max.year.sampling,
                       sampling.years = sampling.years,
                       N.sampling.years = N.sampling.years,
                       age.min = age.min,
                       age.max = age.max,
                       N.FOI = N.FOI,
                       N.birth.years = N.birth.years,
                       birth.years = birth.years,
                       N.titer.sets = N.titer.sets,
                       titer.observable.max = titer.observable.max,
                       N.titers = N.titers,
                       Titers.0 = Titers.0)



# Analysis of seroprevalence data to compare with the model with antibody titers
library(rstan)
## plot seroprevalence data ----
data_1 = get_data_seroprevalence(data.EV71.Malaysia,titer.threshold = 1)
data_2 = get_data_seroprevalence(data.EV71.Malaysia,titer.threshold = 2)
data_3 = get_data_seroprevalence(data.EV71.Malaysia,titer.threshold = 3)

plot_seroprevalence_data(data_1)
dev.copy(pdf,"results/Seroprevalence_1.pdf", width = 7, height = 5)
dev.off()

plot_seroprevalence_data(data_2)
dev.copy(pdf,"results/Seroprevalence_2.pdf", width = 7, height = 5)
dev.off()

plot_seroprevalence_data(data_3)
dev.copy(pdf,"results/Seroprevalence_3.pdf", width = 7, height = 5)
dev.off()



titer.threshold = 3
if(titer.threshold == 1){
  data = data_1
}
if(titer.threshold == 2){
  data = data_2
}
if(titer.threshold == 3){
  data = data_3
}
n.samples <- seropositive <- matrix(data = 0, nrow = N.sampling.years, ncol = age.max)
for(s in seq(min.year.sampling, max.year.sampling) ){
  for(a in seq(age.min, age.max)){
    d = data %>% filter(age == a & sampling.year ==s )
    n.samples[s-min.year.sampling+1, a-age.min+1] = d$n.samples
    seropositive [s-min.year.sampling+1, a-age.min+1] = d$seropositive
  }
}


#  to do: implement in rstan a serocatalytic model to estimate the FOI between 1983 and 2012 and the decay rate omega
data.list = list(NFOI = N.FOI,
                 Nyears = N.sampling.years,
                 NAges = age.max,
                 nsamples = n.samples,
                 seropositive = seropositive)


model = rstan::stan_model(file  = "R/seroprevalence_model.stan")
model.1993 = rstan::stan_model(file  = "R/seroprevalence_model_1993.stan")

rstan.fit = rstan::sampling(object = model, data.list)
rstan.fit.1993 = rstan::sampling(object = model.1993, data.list)

plot_serocatalytic_foi(rstan.fit)

#dev.copy(pdf,paste0("results/Attack_rate_",titer.threshold,".pdf"), width = 5, height = 4)
#dev.off()

plot_serocatalytic_fit(rstan.fit.1993,data)
dev.copy(pdf,paste0("results/Seroprevalence_fit_",titer.threshold,".pdf"), width = 7, height = 5)
dev.off()

saveRDS(object = rstan.fit, file = paste0("results/Serocatalytic_threshold_", titer.threshold,'.rds') )
saveRDS(object = rstan.fit.1993, file = paste0("results/Serocatalytic_1993_threshold_", titer.threshold,'.rds') )


## age at first infection by birth year -----
titer.threshold=3

if(titer.threshold == 1){
  data = data_1
}
if(titer.threshold == 2){
  data = data_2
}
if(titer.threshold == 3){
  data = data_3
}
n.samples <- seropositive <- matrix(data = 0, nrow = N.sampling.years, ncol = age.max)
for(s in seq(min.year.sampling, max.year.sampling) ){
  for(a in seq(age.min, age.max)){
    d = data %>% filter(age == a & sampling.year ==s )
    n.samples[s-min.year.sampling+1, a-age.min+1] = d$n.samples
    seropositive [s-min.year.sampling+1, a-age.min+1] = d$seropositive
  }
}

data.list = list(NFOI = N.FOI,
                 Nyears = N.sampling.years,
                 NAges = age.max,
                 nsamples = n.samples,
                 seropositive = seropositive)


rstan.fit = readRDS(file = paste0("results/Serocatalytic_threshold_", titer.threshold,'.rds') )
rstan.fit.1993 = readRDS(file = paste0("results/Serocatalytic_1993_threshold_", titer.threshold,'.rds') )

fit = rstan.fit.1993
fit = rstan.fit


Chains=rstan::extract(fit)

foi  = 1-exp(-colMeans(Chains$FOI[,1:N.FOI]))
df = NULL
I=0
for(birth in c( 1995,2000,2005)){
  I=I+1

  J = which(birth.years == birth)
  FOI = foi[seq(J,min(J+11, N.FOI))]
  A  = length(FOI)
  C = rep(0,A)
  C[1] = FOI[1]
  if(A>1){
    for(a in 2:A){
      P=1
      for(j in 1:(a-1)){
        P = P*(1-FOI[j])
      }
      C[a] = FOI[a]*P
    }
  }


  df=rbind(df, data.frame(age= 1:A, C = C, S= birth))
}

g= df %>%
  mutate(Year=as.factor(S)) %>%
  ggplot() +
  geom_col(aes(x=age, y = C, group =Year, fill=Year),position = "dodge")+
  theme_bw()+
  ylab('Proportion')+
  xlab('Age at first infection')+
  labs(fill ="Birth year")+
  scale_x_continuous(labels = as.character(1:age.max), breaks = seq(1,age.max))+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)
dev.copy(pdf,paste0("results/Age_infection_serocatalytic_",titer.threshold,".pdf"), width = 6, height = 4)
dev.off()


# Mean age at first infection ----
df = NULL
I=0
for(birth in seq(1995, 2005)){
  I=I+1
  J = which(birth.years == birth)

  FOI = foi[seq(J,min(J+11, N.FOI))]
  A  = length(FOI)
  C = rep(0,A)
  C[1] = FOI[1]
  if(A>1){
    for(a in 2:A){
      P=1
      for(j in 1:(a-1)){
        P = P*(1-FOI[j])
      }
      C[a] = FOI[a]*P
    }
  }
  df=rbind(df, data.frame(age= 1:A, C = C, S= birth))
}
A=df %>%
  mutate(Year=as.factor(S)) %>%
  group_by(Year) %>%
  summarise(a = sum(age*C))


g= df %>%
  group_by(S) %>%
  summarise(a = sum(age*C)) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x=S, y = a-1), size=1.2)+
  theme_bw()+
  ylab('Mean age at first infection')+
  xlab('Birth year')+
  ylim(c(0,5)) +
  scale_x_continuous(labels = as.character(c(1995,2000,2005)), breaks = c(1995,2000,2005))+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)

dev.copy(pdf,paste0("results/Mean_Age_infection_serocatalytic_",titer.threshold,".pdf"), width = 4, height = 4)
dev.off()



 n.sim=100
 foi  = 1-exp(-Chains$FOI[,1:N.FOI])


#foi  = 1-exp(-(res$params[5000:10000,1:N.FOI]))
indices = sample(x = seq(1,nrow(foi)), n.sim)

g= indices%>%
  map(probability_first_infection) %>%
  bind_rows() %>%
  group_by(year_first_infection) %>%
  # summarise(age_mean = sum((age-1)*C)/sum(C)) %>%
  filter(year_first_infection >=1995) %>%
  summarise_at(.vars = "age_mean",
               .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
  ggplot() +
  geom_line(aes(x=year_first_infection, y = mean), size=1.2)+
  geom_ribbon(aes(x = year_first_infection, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= "black") +
  theme_bw()+
  ylab('Mean age at first infection')+
  xlab('Year')+
  scale_x_continuous(labels = as.character(c(1995,2000,2005, 2010)), breaks = c(1995,2000,2005, 2010))+
  ylim(c(0,4)) +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)

dev.copy(pdf,paste0("results/Mean_Age_infection_by_year_serocatalytic_",titer.threshold,".pdf"), width = 4, height = 4)
dev.off()






## decay time ----

decay_curve <- function(omega){
  Time = seq(0,8,by=0.1)
  df  = data.frame(Time = Time, C = exp(-omega*Time))
  return(df)
}


omega = Chains$omega
g = omega  %>%
  map(decay_curve) %>%
  bind_rows() %>%
  group_by(Time) %>%
  summarise_at(.vars = "C",
               .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
  ggplot()+
  geom_line(aes(x=Time, y = mean),col ='dodgerblue3', size=1.2)+
  theme_bw()+
  geom_ribbon(aes(x = Time, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= "dodgerblue3") +
  ylab('Seropositive (%)')+
  xlab('Time (years)')+
  ylim(c(0,1))+
  scale_x_continuous(labels = as.character(0:8), breaks = seq(0,8))+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)
dev.copy(pdf,paste0("results/Survival_titer_",titer.threshold,".pdf"), width = 4, height = 4)
dev.off()


mean(log(2)/omega)

## DIC -----

for(f in c(1,2)){

  if(f==1){
    fit = rstan.fit
  }
  if(f==2){
    fit = rstan.fit.1993
  }

  Chains=rstan::extract(fit)

  P = Chains$P
  dim(Chains$P)
  LL=rep(0,4000)
  Nyears=18
  NAges=12
  for(I in 1:4000){

    ll = 0
    for(j in 1:Nyears){
      for(i in 1:NAges){
        ll = ll+ log(dbinom(x = data.list$seropositive[j,i], size = data.list$nsamples[j,i], prob =  P[I,j,i]))
      }
    }
    LL[I] = ll
  }
  avg_deviance = mean(-2*LL)

  ll.mean = 0
  for(j in 1:Nyears){
    for(i in 1:NAges){
      ll.mean = ll.mean+ log(dbinom(x = data.list$seropositive[j,i], size = data.list$nsamples[j,i], prob =  mean(P[,j,i])))
    }
  }
  deviance_avg_params = -2*ll.mean

  pD = avg_deviance-deviance_avg_params
  DIC = pD + avg_deviance
  print(DIC)

}


