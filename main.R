rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)

source("R/Model_EV71.R")
source("R/utils.R")


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
inds_to_update <- 1:length(params0) # Here we update all the parameters
inds_to_update <- c(31,32) # Here we update all the parameters

mcmc_steps <- 70000
mcmc_adaptive_steps <- 20000


mcmc_steps <- 50
mcmc_adaptive_steps <- 40

is_invalid <- function(k, value) { # Function that checks if parameter value is invalid
  if (k == (n_change_points + 2) & value > 1) { return(TRUE) } # Prop VOC at t=0 cannot be > 1
  FALSE
}

is_invalid <- function(k, value) { # Function that checks if parameter value is invalid
  if (value <10-10) { return(TRUE) } # All the parameters must be > 0
  if (k == 31 & value >8) { return(TRUE) }  
  if (k == 32 & value >8) { return(TRUE) }  
  FALSE
}


res <- run_MCMC(compute_loglik, is_invalid, params0,
                inds_to_update = inds_to_update,
                mcmc_steps = mcmc_steps,
                mcmc_adaptive_steps = mcmc_adaptive_steps,
                verbose = TRUE)

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
