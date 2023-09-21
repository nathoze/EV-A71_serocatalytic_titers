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






data = data_3
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

rstan_fit = rstan::sampling(object = model, data.list)

plot_serocatalytic_foi(rstan_fit)

dev.copy(pdf,"results/Attack_rate_3.pdf", width = 8, height = 4)
dev.off()

plot_serocatalytic_fit(rstan_fit,data)

dev.copy(pdf,"results/Seroprevalence_fit_3.pdf", width = 7, height = 5)
dev.off()

