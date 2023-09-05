# Analysis of seroprevalence data to compare with the model with antibody titers
library(rstan)
## plot seroprevalence data ----
data_1 = get_data_seroprevalence(data.EV71.Malaysia,titer.threshold = 1)
data_2 = get_data_seroprevalence(data.EV71.Malaysia,titer.threshold = 2)
data_3 = get_data_seroprevalence(data.EV71.Malaysia,titer.threshold = 3)
plot_seroprevalence_data(data_1)
plot_seroprevalence_data(data_2)
plot_seroprevalence_data(data_3)

data = data_1
n.samples <- seropositive <- matrix(data = 0, nrow = N.sampling.years, ncol = age.max)
for(s in seq(min.year.sampling, max.year.sampling) ){
  for(a in seq(age.min, age.max)){
    d = data %>% filter(age == a & sampling.year ==s )
    n.samples[s-min.year.sampling+1, a-age.min+1] = d$n.samples
    seropositive [s-min.year.sampling+1, a-age.min+1] = d$seropositive
  }
}


#  to do: implement in rstan a serocatalytic model to estimate the FOI between 1983 and 2012 and the decay rate omega
data= list(NFOI = N.FOI,
           Nyears = N.sampling.years,
           NAges = age.max,
           nsamples = n.samples,
           seropositive = seropositive)


model = rstan::stan_model(file  = "R/seroprevalence_model.stan")

fit = rstan::sampling(object = model,data)
