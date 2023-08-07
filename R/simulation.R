
d1 = data.EV71.Malaysia %>%
  group_by(age, sampling.year) %>%
  summarise(N = sum(n))

# generate data with the same distribution for the total number of samples for a given
# age class and sampling year, with N>0
mu = mean(d1$N)-1
p = mu/var(d1$N)
n=nrow(d1)
RNB = rnbinom(n = nrow(d1), mu = mu, size=mu*p/(1-p)) +1



results = readRDS(file='results/Model_constant.rds')
params = results$params[nrow(results$params),]

all.params = results$model$get_all_parameters(params = params,
                                              fct_model_antibody_increase = results$model$fct_model_antibody_increase,
                                              fct_model_antibody_decrease = results$model$fct_model_antibody_decrease)
data_simulated = all.params$titer.distribution %>%
  mutate(N = rep(RNB, each=5)) %>%
  group_by(age, sampling.year) %>%
  mutate(n = rmultinom(n = 1, size = N,prob = obs.proportion )) %>%
  select(sampling.year, titer.class, birth.year,  n, age)



min.year.sampling <- min(data_simulated$sampling.year)
max.year.sampling <- max(data_simulated$sampling.year)
sampling.years = seq(min.year.sampling, max.year.sampling)
N.sampling.years = length(sampling.years)
age.min = 1 # minimal age in all the surveys at the time of survey
age.max = 12 # maximal age in all surveys at the time of survey
N.FOI = N.sampling.years + age.max # The number of FOI measured. We consider here that it is possible to be infected at age 0
N.birth.years   = N.FOI # The number of birth years considered (-1 if we remove the latest year)
birth.years = seq(min.year.sampling-age.max, max.year.sampling)

#start.sampling.year = age.max
N.titer.sets = age.max*N.sampling.years
titer.observable.max = max(data_simulated$titer.class)

N.titers = 10 # the number of possible values of the titers
Titers.0 <- c(1,rep(0,N.titers-1)) # probability distribution of the titers at birth (everybody is seronegative at birth)




################################################################################
# MCMC
################################################################################

mcmc_steps <- 10000
mcmc_adaptive_steps <- 5000

# Example 1: Constant Model for the FOI ----

params0=c(0.4427671, 1.1189196, 1.73814652)
inds_to_update <- 1:length(params0)
model_constant  =  define_model(fct_model_antibody_increase = get_increase_matrix,
                                fct_model_antibody_decrease = get_decay_matrix,
                                compute_loglik = compute_loglik,
                                params0 = params0,
                                inds_to_update = inds_to_update,
                                model_foi = 'constant')


res <-  run_MCMC_specify_model(model = model_constant,
                               data = data_simulated,
                               mcmc_steps = mcmc_steps,
                               mcmc_adaptive_steps = mcmc_adaptive_steps,
                               verbose = TRUE)

saveRDS(res, file='results/Simulated_data_model_constant.rds')


plot_fit(res, burn_in = 50,n.sim = 10)
