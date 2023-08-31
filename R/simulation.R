min.year.sampling <- 1995
max.year.sampling <- 2012
sampling.years = seq(min.year.sampling, max.year.sampling)
N.sampling.years = length(sampling.years)
age.min = 1 # minimal age in all the surveys at the time of survey
age.max = 12 # maximal age in all surveys at the time of survey
N.FOI = N.sampling.years + age.max # The number of FOI measured. We consider here that it is possible to be infected at age 0
N.birth.years   = N.FOI # The number of birth years considered (-1 if we remove the latest year)
birth.years = seq(min.year.sampling-age.max, max.year.sampling)

N.titer.sets = age.max*N.sampling.years
titer.observable.max = 5

N.titers = 10 # the number of possible values of the titers
Titers.0 <- c(1,rep(0,N.titers-1)) # probability distribution of the titers at birth (everybody is seronegative at birth)


# ## FOI model for the simulated dataset
# model.constant  =  define_model(fct_model_antibody_increase = get_increase_matrix,
#                                 fct_model_antibody_decrease = get_decay_matrix,
#                                 compute_loglik = compute_loglik,
#                                 params0 = params0,
#                                 inds_to_update = inds_to_update,
#                                 model_foi = 'constant')
# params= c(0.45,1.2,1.7)
#
# #results = readRDS(file='results/Model_constant.rds')
# #params = results$params[nrow(results$params),]
#
# all.params = model.constant$get_all_parameters(params = params,
#                                                fct_model_antibody_increase = model.constant$fct_model_antibody_increase,
#                                                fct_model_antibody_decrease = model.constant$fct_model_antibody_decrease)
#
#
# # generate data with the same distribution for the total number of samples for a given
# # age class and sampling year, with N>0
# d1 = data.EV71.Malaysia %>%
#   group_by(age, sampling.year) %>%
#   summarise(N = sum(n))
#
# mu = mean(d1$N)-1
# p = mu/var(d1$N)
# n=nrow(d1)
# RNB = rnbinom(n = nrow(d1), mu = mu, size=mu*p/(1-p)) +1
#
#
# data.simulated = all.params$titer.distribution %>%
#   mutate(N = rep(RNB, each = titer.observable.max)) %>%
#   group_by(age, sampling.year) %>%
#   mutate(n = rmultinom(n = 1, size = N,prob = obs.proportion )) %>%
#   select(sampling.year, titer.class, birth.year,  n, age)

simulate_dataset <- function(model){
  all.params = model$get_all_parameters(params = model$params0,
                                        fct_model_antibody_increase = model$fct_model_antibody_increase,
                                        fct_model_antibody_decrease = model$fct_model_antibody_decrease)

  # generate data with the same distribution for the total number of samples for a given
  # age class and sampling year, with N>0
  d1 = data.EV71.Malaysia %>%
    group_by(age, sampling.year) %>%
    summarise(N = sum(n))

  mu = mean(d1$N)#-1
  p = (mu)/var(d1$N)
  n=nrow(d1)
  RNB = rnbinom(n = nrow(d1), mu = mu, size=(mu+1)*p/(1-p)) #+1


  data.simulated = all.params$titer.distribution %>%
    mutate(N = rep(RNB, each = titer.observable.max)) %>%
    group_by(age, sampling.year) %>%
    mutate(n = rmultinom(n = 1, size = N,prob = obs.proportion )) %>%
    select(sampling.year, titer.class, birth.year,  n, age)

  return(data.simulated)
}


params0=c(0.05,20,3,1.1,2)
n_params = length(params0)
inds_to_update <- 1:length(params0) # Here we update all the parameters
model.peak.constant =  define_model(fct_model_antibody_increase = get_increase_matrix,
                                    fct_model_antibody_decrease = get_decay_matrix,
                                    compute_loglik = compute_loglik,
                                    params0 = params0,
                                    inds_to_update = inds_to_update,
                                    model_foi = "peak_constant")

data.simulated <- simulate_dataset(model.peak.constant)

plot_data(data.simulated)


## MCMC ----

mcmc_steps <- 10000
mcmc_adaptive_steps <- 5000

# Example 1: Constant Model for the FOI ----

params0=c(0.4427671, 1.1189196, 1.73814652)
inds_to_update <- 1:length(params0)
model.constant  =  define_model(fct_model_antibody_increase = get_increase_matrix,
                                fct_model_antibody_decrease = get_decay_matrix,
                                compute_loglik = compute_loglik,
                                params0 = params0,
                                inds_to_update = inds_to_update,
                                model_foi = 'constant')


res <-  run_MCMC(model = model.constant,
                 data = data.simulated,
                 mcmc_steps = mcmc_steps,
                 mcmc_adaptive_steps = mcmc_adaptive_steps,
                 verbose = TRUE)

saveRDS(res, file='results/Simulated_data_model_constant.rds')


## Analyse results of the inference -----
res = readRDS(file = 'results/Simulated_data_model_constant.rds')
data.simulated= res$data
plot_fit(res, burn_in = 6000,n.sim = 30)

hist(res$params[,1])
hist(res$params[,2])
hist(res$params[,3])
