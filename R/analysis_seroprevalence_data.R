# Analysis of seroprevalence data to compare with the model with antibody titers

## plot seroprevalence data ----
data_1 = get_data_seroprevalence(data.EV71.Malaysia,titer.threshold = 1)
data_2 = get_data_seroprevalence(data.EV71.Malaysia,titer.threshold = 2)
data_3 = get_data_seroprevalence(data.EV71.Malaysia,titer.threshold = 3)
plot_seroprevalence_data(data_1)
plot_seroprevalence_data(data_2)
plot_seroprevalence_data(data_3)

#  to do: implement in rstan a serocatalytic model to estimate the FOI between 1983 and 2012 and the decay rate omega
data= list(N.FOI = N.FOI,
           Nyears = N.sampling.years,
           NAge = age.max)
