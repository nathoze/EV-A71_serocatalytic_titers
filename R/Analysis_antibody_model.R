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
#N.FOI = N.sampling.years + age.max # The number of FOI measured. We consider here that it is possible to be infected at age 0
N.FOI = N.sampling.years + age.max-age.min # The number of FOI measured. We consider here that it is possible to be infected at age 0
N.birth.years   = N.FOI # The number of birth years considered (-1 if we remove the latest year)
birth.years = seq(min.year.sampling-age.max, max.year.sampling)

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

# Model comparison -----

C =rep(0,15)
res= readRDS(file='results/Model_constant_no_protection.rds')
C[1] = compute_DIC(res, burn_in = 5000)$DIC
res = readRDS(file='results/Model_5years_no_protection.rds')
C[2] = compute_DIC(res, burn_in = 5000)$DIC
res= readRDS(file='results/Model_independent_no_protection.rds')
C[3] = compute_DIC(res, burn_in = 5000)$DIC
res= readRDS(file='results/Model_peak_constant_no_protection.rds')
C[4] = compute_DIC(res, burn_in = 5000)$DIC
res= readRDS(file='results/Model_independent_1993_no_protection.rds')
C[5] = compute_DIC(res, burn_in = 5000)$DIC

res= readRDS(file='results/Model_constant_no_seroreversion.rds')
C[6] = compute_DIC(res, burn_in = 5000)$DIC
res = readRDS(file='results/Model_5years_no_seroreversion.rds')
C[7] = compute_DIC(res, burn_in = 5000)$DIC
res= readRDS(file='results/Model_independent_no_seroreversion.rds')
C[8] = compute_DIC(res, burn_in = 5000)$DIC
res= readRDS(file='results/Model_peak_constant_no_seroreversion.rds')
C[9] = compute_DIC(res, burn_in = 5000)$DIC
res= readRDS(file='results/Model_independent_1993_no_seroreversion.rds')
C[10] = compute_DIC(res, burn_in = 5000)$DIC

res= readRDS(file='results/Model_constant_protection.rds')
C[11] = compute_DIC(res, burn_in = 5000)$DIC
res = readRDS(file='results/Model_5years_protection.rds')
C[12] = compute_DIC(res, burn_in = 5000)$DIC
res= readRDS(file='results/Model_independent_protection.rds')
C[13] = compute_DIC(res, burn_in = 5000)$DIC
res= readRDS(file='results/Model_peak_constant_protection.rds')
C[14] = compute_DIC(res, burn_in = 5000)$DIC
res= readRDS(file='results/Model_independent_1993_protection.rds')
C[15] = compute_DIC(res, burn_in = 5000)$DIC


plot_fit(res, burn_in = 5000, n.sim = 30, individual.points = TRUE)
plot_foi(res,burn_in = 5000, n.sim=200, show.attack.rate = TRUE)
plot_foi(res,burn_in = 5000, n.sim=200, show.attack.rate = FALSE)
#
# mcmc_steps=10000
# burn_in <- 5000
# thinning <- seq(burn_in + 1, mcmc_steps, by = 5)
# chain <- res$params[thinning, ]
#
# # Visualize MCMC chain
# #colnames(chain) <- c("R0 hist", "pVOC_0", "Advantage")
# chain %>%
#   as_tibble() %>%
#   mutate(iter = seq_len(n())) %>%
#   gather(-iter, key = "variable", value = "value") %>%
#   ggplot() +
#   geom_line(aes(x = iter, y = value, color = variable)) +
#   facet_wrap(~variable, scale = "free") +
#   theme_bw() + theme(legend.position = "none")
#
# res= readRDS(file='results/Model_constant_protection.rds')
# compute_DIC(res, burn_in = 5000)
#
#
# res= readRDS(file='results/Model_constant_no_protection.rds')
# compute_DIC(res, burn_in = 5000)
#

## compare results ----


# res = readRDS(file='results/Model_5years_no_seroreversion.rds')
# print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
# res= readRDS(file='results/Model_independent_no_seroreversion.rds')
# print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
# res= readRDS(file='results/Model_constant_no_seroreversion.rds')
# print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
# res= readRDS(file='results/Model_peak_constant_no_seroreversion.rds')
# print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
#
# res = readRDS(file='results/Model_5years_protection.rds')
# print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
# res= readRDS(file='results/Model_independent_protection.rds')
# print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
# res= readRDS(file='results/Model_constant_protection.rds')
# print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
# res= readRDS(file='results/Model_peak_constant_protection.rds')
# print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )


#  The two best models are the independent models ----
res= readRDS(file='results/Model_independent_protection.rds')
compute_DIC(res, burn_in = 5000)
alpha = res$params[,N.FOI+3]
P = 100*(1-exp(-alpha))
100*mean(P)
quantile(100*P,probs = c(0.025,0.975))
print(paste0( round(quantile025(P), digits=3),'  ', round(mean(P), digits=3), '  ',round(quantile975(P) , digits=3)) )



res= readRDS(file='results/Model_independent_no_protection.rds')
compute_DIC(res, burn_in = 5000)

plot_foi(res,burn_in = 5000, n.sim=200, show.attack.rate = TRUE)
dev.copy(pdf,"results/Attack_Rate_antibody_model.pdf", width = 8, height = 4)
dev.off()


# plot results : the response of a naive individual ----


Titer.max = 5

lambda = res$params[5000:10000,N.FOI+1]
plot(cumsum(truncated_poisson_pmf(lambda = lambda[1], truncation.point = N.titers ,1:11)))
# Average titer increase :
print(paste0( round(quantile025(lambda), digits=3),'  ', round(mean(lambda), digits=3), '  ',round(quantile975(lambda) , digits=3)) )


cumulative_response_distribution <- function(index){
  C = c(0, cumsum(truncated_poisson_pmf(lambda = lambda[index], truncation.point = N.titers ,1:(Titer.max-1))))
  return(data.frame(response = C, titer = 0:(Titer.max-1), titer.exp = 4*2^(0:(Titer.max-1))))
}

response_distribution <- function(index){
  C = c(0, truncated_poisson_pmf(lambda = lambda[index], truncation.point = N.titers ,1:(Titer.max-1)))
  return(data.frame(response = C, titer = 0:(Titer.max-1), titer.exp = 4*2^(0:(Titer.max-1))))

}

n.sim=50
burn_in = 5000
lambda = res$params[-seq(1,burn_in),N.FOI+1]
indices = sample(x = seq(1,length(lambda)), n.sim)

g = indices %>%
  map(response_distribution) %>%
  bind_rows() %>%
  group_by(titer)%>%
  summarise_at(.vars = "response",
               .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
  mutate(titer.exp = 4*2^titer) %>%
  ggplot()+
  # geom_ribbon(aes(x=titer, ymin=quantile025, ymax=quantile975),fill="red", alpha=0.2) +
  #geom_line(aes(x= titer, y= mean), color='red')+
  geom_bar(aes(x= titer, y= mean),stat = "identity", fill  = 'brown2', color='black')+
  xlab('Titer')+
  ylim(c(0, 1))+
  scale_x_continuous(labels = as.character(4*2^seq(0,(Titer.max-1))), breaks = seq(0,(Titer.max-1)))+
  theme_bw()+
  ylab('Frequency')+
  ggtitle('Response after a first infection')+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)
dev.copy(pdf,"results/distribution_response.pdf", width = 3, height = 4)
dev.off()

# plot decay of an individual at 256 (= 7)

decay_distribution <- function(index){
  A = get_decay_matrix(N.titers = 10, omega = omega[index])[Titer.max,]
  C = A[1:(Titer.max)]
  return(data.frame(response = C, titer = 0:(Titer.max-1), titer.exp = 4*2^(0:(Titer.max-1))))
}


omega = res$params[-seq(1,burn_in),N.FOI+2]
indices = sample(x = seq(1,length(omega)), n.sim)

g = indices %>%
  map(decay_distribution) %>%
  bind_rows() %>%
  group_by(titer)%>%
  summarise_at(.vars = "response",
               .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
  mutate(titer.exp = 4*2^titer) %>%
  ggplot()+
  # geom_ribbon(aes(x=titer, ymin=quantile025, ymax=quantile975),fill="red", alpha=0.2) +
  #geom_line(aes(x= titer, y= mean), color='red')+
  geom_bar(aes(x= titer, y= mean),stat = "identity", fill  = 'dodgerblue3', color='black')+
  xlab('Titer')+
  ylim(c(0, 1))+
  scale_x_continuous(labels = as.character(4*2^seq(0,(Titer.max-1))), breaks = seq(0,(Titer.max-1)))+
  theme_bw()+
  ylab('Frequency')+
  ggtitle('Decay')+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)
dev.copy(pdf,"results/distribution_decay.pdf", width = 3, height = 4)
dev.off()

# In the absence of infection, the titer remain stable during the year in xx% of the cases :
d=  indices %>%
  map(decay_distribution) %>%
  bind_rows() %>%
  group_by(titer) %>%
  filter(titer ==4) %>%
  summarise(mean = mean(response), p1 = quantile(response, probs = 0.025), p2 = quantile(response, probs = 0.975))
print(d)

## Cohorts -----

source('R/plot_cohorts.R')


## Comparison results and data by age group ----

plot_titer_age_group_simulation(res)
dev.copy(pdf,"results/Titers_age_class.pdf", width = 14, height = 5)
dev.off()

plot_seroprevalence_age_group_simulation(res)

dev.copy(pdf,"results/Seroprevalence_age_class.pdf", width = 14, height = 5)
dev.off()



## Plot the number of past infections, by age and sampling year  ----

foi  = 1-exp(-colMeans(res$params[5000:10000,1:N.FOI]))
df = NULL

for(sampling.year in c(2000, 2005, 2010)){

  w = which(sampling.year == birth.years)
  foi.years = foi[seq(w-1, w-age.max, by =-1)]
  df=rbind(df, data.frame(age= 1:age.max, C = cumsum(foi.years), S= sampling.year))
}

g= df %>% mutate(Year=as.factor(S)) %>%
  ggplot() + geom_line(aes(x=age, y = C,group =Year, color=Year), size =1.4)+
  theme_bw()+
  ylab('Mean cumulative number of infections')+
  xlab('Age')+
  scale_x_continuous(labels = as.character(1:age.max), breaks = seq(1,age.max))+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)
dev.copy(pdf,"results/Number_infection.pdf", width = 4, height = 4)
dev.off()
print(g)

df = NULL
for(sampling.year in seq(1995, 2011)){

  w = which(sampling.year == birth.years)
  foi.years = foi[seq(w-1, w-age.max, by =-1)]
  df=rbind(df, data.frame(age= 1:age.max, C = cumsum(foi.years), S= sampling.year))
}

g= df %>% mutate(Year=as.factor(S)) %>%
  ggplot() + geom_line(aes(x=age, y = C,group =Year, color=Year), size =1.2)+
  theme_bw()+
  ylab('Mean cumulative number of infections')+
  xlab('Age')+
  scale_x_continuous(labels = as.character(1:age.max), breaks = seq(1,age.max))+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)
dev.copy(pdf,"results/Number_infection_2.pdf", width = 7.5, height = 5)
dev.off()

## Distribution of the age at first infection  ----

foi  = 1-exp(-colMeans(res$params[5000:10000,1:N.FOI]))
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
dev.copy(pdf,"results/Age_infection.pdf", width = 6, height = 4)
dev.off()



# Mean age at first infection by birth year ----

foi  = 1-exp(-colMeans(res$params[5000:10000,1:N.FOI]))
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

g= df %>%
  group_by(S) %>%
  summarise(a = sum(age*C)) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x=S, y = a-1), size=1.2)+
  theme_bw()+
  ylab('Mean age at first infection')+
  xlab('Birth year')+
  scale_x_continuous(labels = as.character(c(1995,2000,2005)), breaks = c(1995,2000,2005))+
  ylim(c(0,NA)) +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)

dev.copy(pdf,"results/Mean_Age_infection.pdf", width = 4, height = 4)
dev.off()




# Mean age at first infection by  year ----
# we have to assume a constant demography
probability_first_infection <- function(index.foi){

  prob =  foi[index.foi,]
  df = NULL
  I=0
  for(birth in seq(1983, 2011)){
    I=I+1
    J = which(birth.years == birth)

    FOI = prob[seq(J,min(J+11, N.FOI))]
    A  = length(FOI)
    age = seq(1,A)
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
    df=rbind(df, data.frame(age= age, C = C, birth= birth, year_first_infection = birth+age-1))
  }

  df= df%>%
    group_by(year_first_infection) %>%
    summarise(age_mean = sum((age-1)*C)/sum(C))

  return(df)
}
n.sim=100
foi  = 1-exp(-(res$params[5000:10000,1:N.FOI]))
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
  #  ylim(c(0,4)) +
  ylim(c(0,2)) +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)
dev.copy(pdf,"results/Mean_Age_infection_by_Year.pdf", width = 5, height = 4)
dev.off()


##  Number of new infections on each year ----
# Number of infections that are new infections

foi  = 1-exp(-(res$params[5000:10000,1:N.FOI]))
foi = foi[1,]
params = res$params
res$model$get_all_parameters(params = params[10,],
                             fct_model_antibody_increase = res$model$fct_model_antibody_increase,
                             fct_model_antibody_decrease = res$model$fct_model_antibody_decrease)

N=100000
population.negative <- rep(N/age.max, age.max)
population.positive <- rep(0, age.max)
sum.new.infections = rep(0, N.FOI) # the first infections
sum.infections = rep(0, N.FOI)
New.Infections = array(data=0, dim=c(N.FOI, age.max))
Infections = array(data=0, dim=c(N.FOI, age.max))
Pop.Neg = array(data=0, dim=c(N.FOI, age.max))
Pop.Pos = array(data=0, dim=c(N.FOI, age.max))

I=0
for(year in birth.years){
  I=I+1
  new.infections = foi[I]*population.negative
  population.positive = population.positive + new.infections
  population.negative = population.negative - new.infections

  population.negative[2:age.max] =  population.negative[seq(1,age.max-1)]
  population.positive[2:age.max] =  population.positive[seq(1,age.max-1)]
  population.negative[1] =N/age.max
  population.positive[1] = 0
  sum.new.infections[I] = sum(new.infections)
  New.Infections[I,] = new.infections
  Pop.Neg[I,] = population.negative
  Pop.Pos[I,] = population.positive
  sum.infections[I] = sum(foi[I]*(population.negative+population.positive))
}

df = data.frame(first.infections =sum.new.infections, all.infections = sum.infections, years = birth.years)

foi  = 1-exp(-(res$params[5000:10000,1:N.FOI]))

incidence_new_infection <- function(index){

  prob.infection = foi[index,]

  N=10000
  population.negative <- rep(N/age.max, age.max)
  population.positive <- rep(0, age.max)
  sum.new.infections = rep(0, N.FOI) # the first infections
  sum.infections = rep(0, N.FOI)
  New.Infections = array(data=0, dim=c(N.FOI, age.max))
  Infections = array(data=0, dim=c(N.FOI, age.max))
  Pop.Neg = array(data=0, dim=c(N.FOI, age.max))
  Pop.Pos = array(data=0, dim=c(N.FOI, age.max))

  I=0
  for(year in birth.years){
    I=I+1
    new.infections = prob.infection[I]*population.negative
    population.positive = population.positive + new.infections
    population.negative = population.negative - new.infections

    population.negative[2:age.max] =  population.negative[seq(1,age.max-1)]
    population.positive[2:age.max] =  population.positive[seq(1,age.max-1)]
    population.negative[1] =N/age.max
    population.positive[1] = 0
    sum.new.infections[I] = sum(new.infections)
    New.Infections[I,] = new.infections
    Pop.Neg[I,] = population.negative
    Pop.Pos[I,] = population.positive
    sum.infections[I] = sum(prob.infection[I]*(population.negative+population.positive))
  }

  df = data.frame(first.infections =sum.new.infections, all.infections = sum.infections, years = birth.years) %>%
    mutate(incidence.new = first.infections/N*10000) %>%
    mutate(incidence.all = all.infections/N*10000)
  return(df)
}

A=indices %>%
  map(incidence_new_infection) %>%
  bind_rows()%>%
  group_by(years) %>%
  filter(years >=1995) %>%
  summarise_at(.vars = c("incidence.new", "incidence.all"),
               .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
  ggplot() +
  geom_line(aes(x=years, y = incidence.new_mean), size=1.2)+
  geom_ribbon(aes(x = years, ymin= incidence.new_quantile025, ymax = incidence.new_quantile975), alpha=0.2, fill= "black") +
  geom_line(aes(x=years, y = incidence.all_mean), color='red', size=1.2)+
  geom_ribbon(aes(x = years, ymin = incidence.all_quantile025, ymax=incidence.all_quantile975), alpha=0.2, fill= "red") +
  #geom_line(aes(x= years, y = incidence.all)) +
  theme_bw()+
  ylab('Incidence per 10,000 children per year')+
  xlab('Year')+
  scale_x_continuous(labels = as.character(c(1995,2000,2005, 2010)), breaks = c(1995,2000,2005, 2010))+
  #  ylim(c(0,4)) +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(A)
dev.copy(pdf,"results/Incidence.pdf", width =5, height = 4)
dev.off()

g = df %>% filter(years >=1995) %>%
  mutate(incidence.new = first.infections/N*10000) %>%
  mutate(incidence.all = all.infections/N*10000) %>%
  ggplot()+
  geom_line(aes(x= years, y = incidence.new))+
  geom_line(aes(x= years, y = incidence.all)) +
  theme_bw()+
  ylab('Incidence per 10,000 children per year')+
  xlab('Year')+
  scale_x_continuous(labels = as.character(c(1995,2000,2005, 2010)), breaks = c(1995,2000,2005, 2010))+
  #  ylim(c(0,4)) +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)






rownames(Pop.Neg) = birth.years
rownames(Pop.Pos) = birth.years
plot(sum.new.infections[13:29]/N*10000,type='l')
plot((New.Infections[13:29,1]),type='l') # infections in young kids
plot((New.Infections[13:29,3]),type='l') # infections in young kids







## Time to reach the minimal titer value ----
n.sim=50
omega = res$params[5000:10000,N.FOI+2]
indices = sample(x = seq(1,length(omega)), n.sim)

survival_64 <- function(omega){
  Mdecay = get_decay_matrix(omega = omega)

  V0 = rep(0,N.titers)
  V0[5] =1
  VN = V0
  C=rep(0,11)
  df = NULL
  j=0
  for(time in 0:N.titers){
    C[time+1] = VN[1]

    VN = t(Mdecay) %*%VN
  }
  df = data.frame(years = 0:N.titers, C = 100*(1-C))
  return( df)
}


g = omega[indices] %>%
  map(survival_64) %>%
  bind_rows() %>%
  group_by(years) %>%
  summarise_at(.vars = "C",
               .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
  filter(years<=8) %>%
  ggplot()+
  geom_line(aes(x=years, y = mean),col ='dodgerblue3', size=1.2)+
  theme_bw()+
  geom_ribbon(aes(x = years, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= "dodgerblue3") +
  ylab('Seropositive (%)')+
  xlab('Time (years)')+
  scale_x_continuous(labels = as.character(0:8), breaks = seq(0,8))+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)
dev.copy(pdf,"results/Survival.pdf", width = 4, height = 4)
dev.off()
