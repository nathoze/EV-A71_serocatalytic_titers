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
#birth.years = seq(min.year.sampling-age.max, max.year.sampling)
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

################################################################################
# MCMC
################################################################################

mcmc_steps <- 120
mcmc_adaptive_steps <- 50

# Example 1: Constant Model for the FOI ----

params0=c(0.4427671, 1.1189196, 1.73814652,0)
inds_to_update <- 1:length(params0)
model.constant  =  define_model(fct_model_antibody_increase = get_increase_matrix,
                                fct_model_antibody_decrease = get_decay_matrix,
                                compute_loglik = compute_loglik,
                                params0 = params0,
                                inds_to_update = inds_to_update,
                                model_foi = 'constant')


res <-  run_MCMC(model = model.constant,
                 data = data.EV71.Malaysia,
                 mcmc_steps = mcmc_steps,
                 mcmc_adaptive_steps = mcmc_adaptive_steps,
                 verbose = TRUE)

saveRDS(res, file='results/Model_constant.rds')

# Example 2: Model with FOI constant in periods of 5 years ----

params0=c(1.2944222, 0.3230752, 0.6266874, 0.4429604, 0.6127554, 0.7394369, 0.9057352, 2.7328078,0)
n_params = length(params0)
inds_to_update <- 1:length(params0) # Here we update all the parameters
model.five.years =  define_model(fct_model_antibody_increase = get_increase_matrix,
                                 fct_model_antibody_decrease = get_decay_matrix,
                                 compute_loglik = compute_loglik,
                                 params0 = params0,
                                 inds_to_update = inds_to_update,
                                 model_foi = "5years")

res <-  run_MCMC(model = model.five.years,
                 data = data.EV71.Malaysia,
                 mcmc_steps = mcmc_steps,
                 mcmc_adaptive_steps = mcmc_adaptive_steps,
                 verbose = TRUE)

saveRDS(res, file='results/Model_5years.rds')


# Example 3: Independent model ----

params0  = c( runif(n  = N.FOI, max = 0.7) , 3,1,0)
inds_to_update <- 1:length(params0)
model.independent =  define_model(fct_model_antibody_increase = get_increase_matrix,
                                  fct_model_antibody_decrease = get_decay_matrix,
                                  compute_loglik = compute_loglik,
                                  params0 = params0,
                                  inds_to_update = inds_to_update,
                                  model_foi = "independent")


res <-  run_MCMC(model = model.independent,
                 data = data.EV71.Malaysia,
                 mcmc_steps = mcmc_steps,
                 mcmc_adaptive_steps = mcmc_adaptive_steps,
                 verbose = TRUE)

saveRDS(res, file='results/Model_independent.rds')


# Example 4: Peak+constant model ----

params0  = c( 0.1,20,0.5, 3,1,0.0)
inds_to_update <- 1:length(params0)
model.peak.constant =  define_model(fct_model_antibody_increase = get_increase_matrix,
                                    fct_model_antibody_decrease = get_decay_matrix,
                                    compute_loglik = compute_loglik,
                                    params0 = params0,
                                    inds_to_update = inds_to_update,
                                    model_foi = "peak_constant")

res <-  run_MCMC(model = model.peak.constant,
                 data = data.EV71.Malaysia,
                 mcmc_steps = mcmc_steps,
                 mcmc_adaptive_steps = mcmc_adaptive_steps,
                 verbose = TRUE)
saveRDS(res, file='results/Model_peak_constant.rds')



# Example 5: Independent model and no circulation before 1993 ----

params0  = c( runif(n  = N.FOI, max = 0.7) , 3,1,1)
inds_to_update <- 1:length(params0)
params0[1:10] = rep(0,10) # no circulation until 1993 (excluded)
inds_to_update=inds_to_update[-seq(1,10)]

model.independent =  define_model(fct_model_antibody_increase = get_increase_matrix,
                                  fct_model_antibody_decrease = get_decay_matrix,
                                  compute_loglik = compute_loglik,
                                  params0 = params0,
                                  inds_to_update = inds_to_update,
                                  model_foi = "independent")

res <-  run_MCMC(model = model.independent,
                 data = data.EV71.Malaysia,
                 mcmc_steps = mcmc_steps,
                 mcmc_adaptive_steps = mcmc_adaptive_steps,
                 verbose = TRUE)

saveRDS(res, file='results/Model_independent_1993_protection.rds')

