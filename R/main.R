rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)

source("R/Model_EV71.R")
source("R/utils.R")
source("R/Additional_models.R")


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
N.FOI = N.sampling.years + age.max # The number of FOI measured. We consider here that it is possible to be infected at age 0
N.birth.years   = N.FOI # The number of birth years considered (-1 if we remove the latest year)
birth.years = seq(min.year.sampling-age.max, max.year.sampling)

#start.sampling.year = age.max
N.titer.sets = age.max*N.sampling.years

titer.observable.max = max(data.EV71.Malaysia$titer.class)

N.titers = 10 # the number of possible values of the titers
Titers.0 <- c(1,rep(0,N.titers-1)) # probability distribution of the titers at birth (everybody is seronegative at birth)



################################################################################
# Log-likelihood
################################################################################

compute_loglik <- function(all.params) {

  prob.titers =  left_join(all.params$titer.distribution, data.EV71.Malaysia,
                           by = c("birth.year", "age", "titer.class", "sampling.year"))  %>%
    group_by(age, sampling.year) %>%
    summarise(ll  = dmultinom(x = n,prob = obs.proportion, log=TRUE), .groups='drop')  %>%
    summarise(s= sum(ll))

  return(prob.titers$s)

}


################################################################################
# MCMC
################################################################################


params0  = c( runif(n  = N.FOI, max = 0.7) , 3,1)
#params0 <- c(c(rep(1.0, n_change_points + 1), 0.02), 0.6)
#n_params <- length(params0)
# # inds_to_update contains indices of parameters to update. Here we do not update
# #   first parameter (R of historical virus), which will be fixed at 1.0
# inds_to_update <- 2:length(params0)
inds_to_update <- c(31,32) # Here we update the antibody model parameters
inds_to_update <- 1:length(params0) # Here we update all the parameters

mcmc_steps <- 120
mcmc_adaptive_steps <- 50

#mcmc_steps <- 50
#mcmc_adaptive_steps <- 40

is_invalid <- function(k, value) { # Function that checks if parameter value is invalid
  if (k == (n_change_points + 2) & value > 1) { return(TRUE) } # Prop VOC at t=0 cannot be > 1
  FALSE
}

is_invalid <- function(k, value) { # Function that checks if parameter value is invalid
  if (value <10e-9) { return(TRUE) } # All the parameters must be > 0
  if (k<=N.FOI & value >2) { return(TRUE) } # the foi
  if (k == N.FOI+1 & value >8) { return(TRUE) }# sigmaP
  if (k == N.FOI+2 & value >8) { return(TRUE) }# Omega
  FALSE
}


res <- run_MCMC(compute_loglik, is_invalid, params0,
                inds_to_update = inds_to_update,
                mcmc_steps = mcmc_steps,
                mcmc_adaptive_steps = mcmc_adaptive_steps,
                verbose = TRUE)


## Define a model
model = list(compute_loglik = compute_loglik,
             is_invalid = is_invalid,
             params0 = params0,
             inds_to_update = inds_to_update,
             get_all_parameters = get_all_parameters)



# Example 1: Specify a  Model constant ----

params0  = c( 0.3, 3,1)
inds_to_update <- 1:length(params0) # Here we update all the parameters

model_constant = list(compute_loglik = compute_loglik,
                      params0 = params0,
                      inds_to_update = inds_to_update,
                      is_invalid = is_invalid_model_constant,
                      get_all_parameters = get_all_parameters_model_constant,
                      update_all_parameters = update_all_parameters_model_constant)

res <-  run_MCMC_specify_model(model = model_constant,
                               mcmc_steps = mcmc_steps,
                               mcmc_adaptive_steps = mcmc_adaptive_steps,
                               verbose = TRUE)

# model_constant$get_all_parameters(res$params[100,])

model_constant$compute_loglik(all.params = model_constant$get_all_parameters(params0))

# Example 2: Constant FOI in periods of five years -----

params0  = c( rep(0.3,6), 3,1)
params0  = c(  runif(n  = round(N.FOI/5), max = 0.7), 3,1)

params0= c( 0.2767452 , 0.2730377, 0.1934079, 0.4318872, 0.3687840, 0.3943001, 1.2745368, 1.0000000)
n_params = length(params0)
inds_to_update <- 1:length(params0) # Here we update all the parameters

model_five_years= list(compute_loglik = compute_loglik,
                       params0 = params0,
                       inds_to_update = inds_to_update,
                       is_invalid = is_invalid_model_five_years,
                       get_all_parameters = get_all_parameters_model_five_years,
                       update_all_parameters = update_all_parameters_model_five_years)

res <-  run_MCMC_specify_model(model = model_five_years,
                               mcmc_steps = mcmc_steps,
                               mcmc_adaptive_steps = mcmc_adaptive_steps,
                               verbose = TRUE)

model_five_years$compute_loglik(all.params = model_five_years$get_all_parameters(params0))


## To do :
# - Need to control the proposal sd because it takes very high values very quickly
# - plot the force of infection per year
# - plot fits of the titer distribution
# - simulate a titer distribution from a model and a set of parameters
# - Implement other models for the titer dynamics and for the FOI


# Plot -----

plot(res$params[1,1:6])

burn_in <- mcmc_adaptive_steps
thinning <- seq(burn_in + 1, mcmc_steps, by = 50)
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

################################################################################
# Simulate trajectories
################################################################################
n_sims <- 500
tf_i_sims <- tf_i_model + 7
chosen_inds <- sample(1:nrow(chain), n_sims)
sims <- matrix(NA, nrow = tf_i_sims, ncol = n_sims)
for (sim in 1:n_sims) {
  chosen_ind <- chosen_inds[sim]
  RHs <- c(chain[chosen_ind, 1])
  pVOC_0 <- chain[chosen_ind, 2]
  alphaVOC <- chain[chosen_ind, 3]
  ps <- pVOC(tf_i_sims, alphaVOC, GT_mu, GT_cv, pVOC_0, RHs, change_points)
  sims[, sim] <- ps
}

sims_plot <- tibble(
  date = t0_model + seq_len(tf_i_sims),
  the_mean = apply(sims, 1, mean),
  lower = apply(sims, 1, quantile, probs = 0.025),
  upper = apply(sims, 1, quantile, probs = 0.975)
)

data_plot <- tibble(
  date = sampling_dates,
  obs = c(data_mat[2] / data_mat[1], data_mat[4] / data_mat[3])
)

sims_plot %>%
  ggplot() +
  geom_ribbon(aes(x = date, ymin = lower, ymax = upper), fill = "dodgerblue3",
              alpha = 0.3) +
  geom_line(aes(x = date, y = the_mean), color = "dodgerblue4") +
  geom_point(data = data_plot, aes(x = date, y = obs), size = 2) +
  xlab("") + ylab("Prop VOC") +
  theme_bw()
