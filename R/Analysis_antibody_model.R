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

# Plot -----
res= readRDS(file='results/Model_constant_no_protection.rds')
compute_DIC(res, burn_in = 5000)
res = readRDS(file='results/Model_5years_no_protection.rds')
compute_DIC(res, burn_in = 5000)
res= readRDS(file='results/Model_independent_no_protection.rds')
compute_DIC(res, burn_in = 5000)
res= readRDS(file='results/Model_peak_constant_no_protection.rds')
compute_DIC(res, burn_in = 5000)
res= readRDS(file='results/Model_independent_1993_no_protection.rds')
compute_DIC(res, burn_in = 5000)

res= readRDS(file='results/Model_constant_protection.rds')
compute_DIC(res, burn_in = 5000)
res = readRDS(file='results/Model_5years_protection.rds')
compute_DIC(res, burn_in = 5000)
res= readRDS(file='results/Model_independent_protection.rds')
compute_DIC(res, burn_in = 5000)
res= readRDS(file='results/Model_peak_constant_protection.rds')
compute_DIC(res, burn_in = 5000)
res= readRDS(file='results/Model_independent_1993_protection.rds')
compute_DIC(res, burn_in = 5000)

res= readRDS(file='results/Model_constant_no_seroreversion.rds')
compute_DIC(res, burn_in = 5000)
res = readRDS(file='results/Model_5years_no_seroreversion.rds')
compute_DIC(res, burn_in = 5000)
res= readRDS(file='results/Model_independent_no_seroreversion.rds')
compute_DIC(res, burn_in = 5000)
res= readRDS(file='results/Model_peak_constant_no_seroreversion.rds')
compute_DIC(res, burn_in = 5000)
res= readRDS(file='results/Model_independent_1993_no_seroreversion.rds')
compute_DIC(res, burn_in = 5000)

plot_fit(res, burn_in = 5000, n.sim = 30, individual.points = TRUE)
plot_foi(res,burn_in = 5000, n.sim=200, show.attack.rate = TRUE)
plot_foi(res,burn_in = 5000, n.sim=200, show.attack.rate = FALSE)

mcmc_steps=10000
burn_in <- 5000
thinning <- seq(burn_in + 1, mcmc_steps, by = 5)
chain <- res$params[thinning, ]

# Visualize MCMC chain
#colnames(chain) <- c("R0 hist", "pVOC_0", "Advantage")
chain %>%
  as_tibble() %>%
  mutate(iter = seq_len(n())) %>%
  gather(-iter, key = "variable", value = "value") %>%
  ggplot() +
  geom_line(aes(x = iter, y = value, color = variable)) +
  facet_wrap(~variable, scale = "free") +
  theme_bw() + theme(legend.position = "none")



res= readRDS(file='results/Model_constant_no_protection.rds')
compute_DIC(res, burn_in = 5000)

res= readRDS(file='results/Model_constant_protection.rds')
compute_DIC(res, burn_in = 5000)

## compare results ----


res = readRDS(file='results/Model_5years_no_seroreversion.rds')
print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
res= readRDS(file='results/Model_independent_no_seroreversion.rds')
print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
res= readRDS(file='results/Model_constant_no_seroreversion.rds')
print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
res= readRDS(file='results/Model_peak_constant_no_seroreversion.rds')
print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )

res = readRDS(file='results/Model_5years_protection.rds')
print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
res= readRDS(file='results/Model_independent_protection.rds')
print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
res= readRDS(file='results/Model_constant_protection.rds')
print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )
res= readRDS(file='results/Model_peak_constant_protection.rds')
print(paste0( round(quantile025(res$params[,ncol(res$params)-1]), digits=2),'  ', round(mean(res$params[,ncol(res$params)-1]), digits=2), '  ',round(quantile975(res$params[,ncol(res$params)-1]) , digits=2)) )


#  Comparative study of the different independent models ----
res= readRDS(file='results/Model_independent_protection.rds')
compute_DIC(res, burn_in = 5000)


res= readRDS(file='results/Model_independent_no_protection.rds')
compute_DIC(res, burn_in = 5000)



plot_foi(res,burn_in = 5000, n.sim=200, show.attack.rate = TRUE)
dev.copy(pdf,"results/Attack_Rate_antibody_model.pdf", width = 8, height = 4)
dev.off()


# plot results : the response of a naive individual ----


Titer.max = 5

lambda = res$params[5000:10000,N.FOI+1]
plot(cumsum(truncated_poisson_pmf(lambda = lambda[1], truncation.point = N.titers ,1:11)))

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
dev.copy(pdf,"results/distribution_response.pdf", width = 4, height = 4)
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
dev.copy(pdf,"results/distribution_decay.pdf", width = 4, height = 4)
dev.off()


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
dev.copy(pdf,"results/Number_infection_2.pdf", width = 6, height = 4)
dev.off()



## Age at first infection ----

foi  = 1-exp(-colMeans(res$params[5000:10000,1:N.FOI]))
df = NULL
foi =0.001*foi


foi[20]=0.5
# problème dans la boucle
for(sampling.year in c(2000, 2005, 2010)){

  w = which(sampling.year == birth.years)
  C = rep(0,age.max)
  C[1] = foi[w-1]
  for(a in 2:age.max){
    P=1
    for(j in 1:(a-1)){
      P = P*(1-foi[w-j])
    }
    C[a] = foi[w-a]*P
  }

   df=rbind(df, data.frame(age= 1:age.max, C = C, S= sampling.year))
}

g= df %>%
  mutate(Year=as.factor(S)) %>%
  ggplot() +
  geom_col(aes(x=age, y = C, group =Year, fill=Year),position = "dodge")+
  theme_bw()+
  ylab('Proportion')+
  xlab('Age at first infection')+
  scale_x_continuous(labels = as.character(1:age.max), breaks = seq(1,age.max))+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        text=element_text(size=18))
print(g)
#dev.copy(pdf,"results/Number_infection_2.pdf", width = 6, height = 4)
#dev.off()
I=0
for(birth in birth.years){
  I=I+1
  FOI = foi[seq(I,min(I+11, N.FOI))]
  A  =length(FOI)
  # C = rep(0,age.max)
  # C[1] = foi[w-1]
  # for(a in 2:age.max){
  #   P=1
  #   for(j in 1:(a-1)){
  #     P = P*(1-foi[w-j])
  #   }
  #   C[a] = foi[w-a]*P
  # }


  # w = which(sampling.year == birth.years)
  # C = rep(0,age.max)
  # C[1] = foi[w-1]
  # for(a in 2:age.max){
  #   P=1
  #   for(j in 1:(a-1)){
  #     P = P*(1-foi[w-j])
  #   }
  #   C[a] = foi[w-a]*P
  # }
  #
  # df=rbind(df, data.frame(age= 1:age.max, C = C, S= sampling.year))
}


