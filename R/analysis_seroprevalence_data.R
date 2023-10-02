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



titer.threshold = 1
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



## DIC serocatalytic model ----

## age at first infection by birth year -----



fit = rstan.fit.1993
fit = rstan.fit
Chains=rstan::extract(fit)


foi  = 1-exp(-colMeans(Chains$FOI[,1:N.FOI]))
df = NULL
I=0
for(birth in c(1990, 1995,2000,2005)){
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



df = NULL
I=0
for(birth in seq(1990, 2005)){
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
  ylim(c(0,NA)) +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)


print(g)
dev.copy(pdf,paste0("results/Mean_Age_infection_serocatalytic_",titer.threshold,".pdf"), width = 4, height = 4)
dev.off()


## Total number of infections by age and by birth year -----



