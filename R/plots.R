# plot the fits of the model to the data
# simulate data from the posterior distribution

simulate_mean_posterior_titer <-function(k, results){

  par <- results$params[k,]
  all.params = results$model$get_all_parameters(params = par,
                                                fct_model_antibody_increase = results$model$fct_model_antibody_increase,
                                                fct_model_antibody_decrease = results$model$fct_model_antibody_decrease)

  td = all.params$titer.distribution

  A = left_join(td, results$data, by = c("birth.year", "age", "titer.class", "sampling.year")) %>%
    group_by(age, sampling.year) %>%
    mutate(N = sum(n)) %>%
    mutate(simul.titer = rmultinom(n = 1, size = N,prob = obs.proportion)*titer.class/N ) %>%
    summarise(mean.titer.sim=colSums(simul.titer)) #%>% mutate(index = k)

  return(A)

}
plot_data <- function(data, individual.points = FALSE){
  g = data%>%
    group_by(age, sampling.year) %>%
    mutate(N = sum(n)) %>%
    mutate(mean.titer.obs = sum(titer.class*n)/N) %>%
    select(sampling.year,age, mean.titer.obs) %>%
    ggplot()+
    geom_line(aes(x=age, y = mean.titer.obs)) +
    geom_point(aes(x=age, y = mean.titer.obs)) +
    facet_wrap(  vars(sampling.year))+
    scale_x_continuous(breaks=seq(1,12))+
    ylab('log2 Titer')+
    ylim(c(0, 5))+
    theme_bw()

  if(individual.points == TRUE){

    data_repeat = data[rep(seq_len(nrow(data)), data$n),]
    g = g +
      geom_jitter(data=data_repeat, aes(x = age , y = titer.class),
                  width = 0.25, height = 0.25, cex=0.3, colour="gray")

  }

  return(g)

}

plot_fit <- function(results, burn_in, n.sim=50, individual.points = FALSE){

  data = results$data%>%
    group_by(age, sampling.year) %>%
    mutate(N = sum(n)) %>%
    mutate(mean.titer.obs = sum(titer.class*n)/N)
  #  %>%   # select(sampling.year,age , mean.titer.obs )

  params = results$params[-seq(1,burn_in),]
  indices = sample(x = seq(1,nrow(params)), n.sim)
  # simulate titers with parameter from the posterior distribution
  g = indices %>%
    map(simulate_mean_posterior_titer, results= results) %>%
    bind_rows() %>%
    left_join(data , by = c( "age", "sampling.year")) %>%
    ggplot()+
    stat_summary(aes(x=age, y = mean.titer.sim),
                 fun.data=mean_sdl, fun.args = list(mult=1),
                 geom="pointrange", color="red")+
    geom_line(aes(x=age, y = mean.titer.obs)) +
    geom_point(aes(x=age, y = mean.titer.obs)) +
    facet_wrap(vars(sampling.year))+
    scale_x_continuous(breaks=seq(1,12))+
    ylim(c(0, 5))+
    ylab('log2 Titer')+
    theme_bw()

  if(individual.points == TRUE){

    data_repeat = data[rep(seq_len(nrow(data)), data$n),]
    g = g +
      geom_jitter(data=data_repeat, aes(x = age , y = titer.class),
                  width = 0.25, height = 0.25, cex=0.3, colour="gray")

  }
  return(g)

}

get_foi <-function(k, results){
  # get the foi by year for a given set of value in the posterior distribution
  params = results$params[k,]
  model_foi = results$model$model_foi

  if(model_foi == 'constant'){
    FOI = rep(params[1], N.FOI)
  }
  if(model_foi == '5years'){
    FOI = rep(0,N.FOI)
    for(i in 1:round(N.FOI/5)){
      FOI[seq((i-1)*5+1,(i)*5)] = params[i]
    }
  }
  if(model_foi == 'independent'){
    FOI = params[1:N.FOI]
  }
  if(model_foi == 'peak_constant'){
    FOI = rep(params[1], N.FOI)
    FOI[round(params[2])] =  FOI[round(params[2])]+params[3]
  }

  df = data.frame(year = birth.years, FOI = FOI)
  return(df)
}

plot_foi <- function(results, burn_in, n.sim = 100, show.attack.rate = TRUE){

  params = results$params
  indices = sample(x = seq(burn_in+1,nrow(params)), n.sim)

  if(show.attack.rate == TRUE){
    ylabel = 'Attack rate'
  }
  if(show.attack.rate == FALSE){
    ylabel = 'Force of infection'
  }
  g = indices %>%
    map(get_foi, results = results) %>%
    bind_rows()%>%
    group_by(year) %>%
    mutate(FOI =  case_when(show.attack.rate ~ 1-exp(-FOI),
                            !show.attack.rate ~ FOI) ) %>%
    summarise_at(.vars = "FOI",
                 .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
    ggplot()+
    geom_line(aes(x= year, y= mean))+
    geom_ribbon(aes(x=year, ymin=quantile025, ymax=quantile975),fill="gray30", alpha=0.2) +
    theme_bw()+
    ylab(ylabel) +
    xlab('')+
    scale_x_continuous(breaks=seq(min(birth.years),max(birth.years), by=2))+
    theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))+
    ylim(c(0,NA))

  return(g)

}

# Plot the fit and the FOI FOI for the simple serocatalytic model ----
# Equivalent to the results above for the titer model
plot_seroprevalence_data <- function(seroprevalence.data){
  g = seroprevalence.data%>%
    group_by(age, sampling.year) %>%
    mutate(lower = binom::binom.confint(seropositive , n.samples, methods = 'exact')$lower) %>%
    mutate(upper = binom::binom.confint(seropositive , n.samples, methods = 'exact')$upper) %>%
    ggplot()+
    geom_line(aes(x=age, y = seroprevalence)) +
    geom_point(aes(x=age, y = seroprevalence)) +
    geom_segment(aes(x = age, y=lower, xend = age, yend = upper)) +
    facet_wrap(  vars(sampling.year))+
    scale_x_continuous(breaks=seq(1,12))+
    ylab('Seroprevalence')+
    ylim(c(0, 1))+
    theme_bw()

  return(g)

}

# simulate_mean_posterior_seroprevalence <- function(k, results){
#
#   FOI <- results$FOI[k,]
#   omega <- results$omega[k]
#   P = matrix(data = 0, nrow = N.sampling.years, ncol = age.max)
#
#   for(birthyear in 1:N.FOI){
#     pinf = 0
#     for(age in 0:age.max){
#       if(age+birthyear <= N.FOI){
#         pinf=pinf*(exp(-FOI[age+birthyear]-omega)) + FOI[age+birthyear]/(FOI[age+birthyear]+omega)*(1-exp(-omega-FOI[age+birthyear]))
#         if(age >0 && (age+birthyear-age.max>= 1)){
#           P[age+birthyear-age.max,age] = pinf
#         }
#       }
#     }
#   }
#   df = data.frame(P  = as.numeric(P),
#                   age = rep(seq(1,age.max), each =N.sampling.years),
#                   sampling.year = rep(sampling.years, age.max) ) %>%
#     group_by(age, sampling.year)
#
#   return(df)
# }




simulate_mean_posterior_seroprevalence <- function(k, results){

  FOI <- results$FOI[k,]
  omega <- results$omega[k]
  P = matrix(data = 0, nrow = N.sampling.years, ncol = age.max)

  for(birthyear in 1:N.FOI){
    pinf = 0
    for(age in 1:age.max){
      if(age+birthyear-1 <= N.FOI){
        pinf=pinf*(exp(-FOI[age+birthyear-1]-omega)) + FOI[age+birthyear-1]/(FOI[age+birthyear-1]+omega)*(1-exp(-omega-FOI[age+birthyear-1]))
        if(age >0 && (age+birthyear-age.max>= 1)){
          P[age+birthyear-age.max,age] = pinf
        }
      }
    }
  }
  df = data.frame(P  = as.numeric(P),
                  age = rep(seq(1,age.max), each =N.sampling.years),
                  sampling.year = rep(sampling.years, age.max) ) %>%
    group_by(age, sampling.year)

  return(df)
}


plot_serocatalytic_fit <- function(rstan_fit, data, n.sim=50){
  Chains = rstan::extract(rstan_fit)
  FOI = Chains$FOI
  omega = Chains$omega
  indices = sample(x = seq(1,nrow(FOI)), n.sim)

  results = list(FOI = FOI,
                 omega = omega)


  g = indices %>%
    map(simulate_mean_posterior_seroprevalence, results = results) %>%
    bind_rows() %>%
    group_by(age,sampling.year) %>%
     summarise_at(.vars = "P",
                 .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
    left_join(data , by = c( "age", "sampling.year")) %>%
    mutate(lower = binom::binom.confint(seropositive , n.samples, methods = 'exact')$lower) %>%
    mutate(upper = binom::binom.confint(seropositive , n.samples, methods = 'exact')$upper) %>%
    ggplot()+
    geom_line(aes(x = age, y = seroprevalence)) +
    geom_point(aes(x = age, y = seroprevalence)) +
    geom_line(aes(x = age, y = mean), color='red')+
    geom_ribbon(aes(x = age, ymin = quantile025, ymax=quantile975), fill="red", alpha=0.2) +
    facet_wrap(vars(sampling.year))+
    scale_x_continuous(breaks=seq(1,12))+
    ylim(c(0, 1))+
    ylab('Seroprevalence')+
    theme_bw()

  return(g)
}

get_serocatalytic_foi <- function(k,results){

  FOI = results[k,]

  return(data.frame(year = birth.years, FOI = FOI))
}

plot_serocatalytic_foi <- function(rstan_fit, n.sim = 100, show.attack.rate = TRUE){

  Chains = rstan::extract(rstan_fit)
  FOI = Chains$FOI

  indices = sample(x = seq(1,nrow(FOI)), n.sim)

  if(show.attack.rate == TRUE){
    ylabel = 'Attack rate'
  }
  if(show.attack.rate == FALSE){
    ylabel = 'Force of infection'
  }

  g = indices %>%
    map(get_serocatalytic_foi, results = FOI) %>%
    bind_rows()%>%
    group_by(year) %>%
    mutate(FOI =  case_when(show.attack.rate ~ 1-exp(-FOI),
                            !show.attack.rate ~ FOI) ) %>%
    summarise_at(.vars = "FOI",
                 .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
    ggplot()+
    geom_line(aes(x= year, y= mean))+
    geom_ribbon(aes(x=year, ymin=quantile025, ymax=quantile975),fill="gray30", alpha=0.2) +
    theme_bw()+
    ylab(ylabel) +
    xlab('')+
    scale_x_continuous(breaks=seq(min(birth.years),max(birth.years), by=2))+
    theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))+
    ylim(c(0,NA))

  return(g)

}



# Data : Mean antibody titers by age group and sampling year ----

plot_titer_age_group <-function(data){

  data.titer <- data%>%
    mutate(age.class =  case_when(age <=3 ~ "1-3",
                                  age >=4 & age <= 6 ~ "4-6",
                                  age> 6 ~ "7-12" ) )%>%
    group_by(age.class, sampling.year) %>%
    mutate(N = sum(n)) %>%
    mutate(mean.titer.obs = sum(titer.class*n)/N) %>%
    select(sampling.year,age.class, N,mean.titer.obs)%>%
    distinct()

  plot.titer.data =  data.titer %>%
    ggplot()+
    geom_line(aes(x= sampling.year, y=mean.titer.obs, group = age.class, color= age.class ))+
    geom_point(aes(x= sampling.year, y=mean.titer.obs, group = age.class, color= age.class, size=N ))+
    theme_bw()+
    ylab("log2 titer") +
    xlab('')+
    ylim(c(1,4.8))+
    scale_x_continuous(breaks=seq(min.year.sampling,max.year.sampling))+
    theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))

  return(plot.titer.data)

}

# Simulation : Mean antibody titers by age group and sampling year ----

simulate_mean_posterior_titer_age_group <-function(k, results){

  par <- results$params[k,]
  all.params = results$model$get_all_parameters(params = par,
                                                fct_model_antibody_increase = results$model$fct_model_antibody_increase,
                                                fct_model_antibody_decrease = results$model$fct_model_antibody_decrease)

  td = all.params$titer.distribution

  A = left_join(td, results$data, by = c("birth.year", "age", "titer.class", "sampling.year")) %>%
    mutate(age.class =  case_when(age <=3 ~ "1-3",
                                  age >=4 & age <= 6 ~ "4-6",
                                  age> 6 ~ "7-12" ) )%>%
    group_by(age.class, sampling.year) %>%
    mutate(N = sum(n)) %>%
    mutate(simul.titer = rmultinom(n = 1, size = N,prob = obs.proportion)*titer.class/N ) %>%
    summarise(mean.titer.sim=colSums(simul.titer)) #%>% mutate(index = k)

  return(A)

}

plot_titer_age_group_simulation <-function(results,burn_in=5000, n.sim = 100){

  data <- results$data

  data.titer <- data%>%
    mutate(age.class =  case_when(age <=3 ~ "1-3",
                                  age >=4 & age <= 6 ~ "4-6",
                                  age> 6 ~ "7-12" ) )%>%
    group_by(age.class, sampling.year) %>%
    mutate(N = sum(n)) %>%
    mutate(mean.titer.obs = sum(titer.class*n)/N) %>%
    select(sampling.year,age.class, N,mean.titer.obs)%>%
    distinct()

  params = results$params[-seq(1,burn_in),]

  indices = sample(x = seq(1,nrow(params)), n.sim)
  d2=indices %>%
    map(simulate_mean_posterior_titer_age_group, results= results) %>%
    bind_rows() %>%
    group_by(age.class, sampling.year) %>%
    summarise_at(.vars = "mean.titer.sim",
                 .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975"))%>%
    left_join(data.titer, by = c( "age.class", "sampling.year") )%>%
    ggplot()+
    geom_line(aes(x= sampling.year, y= mean, group = age.class, color = age.class), linetype='dashed')+
    geom_ribbon(aes(x = sampling.year, ymin=quantile025, ymax=quantile975,group = age.class,fill=age.class), alpha=0.2) +
    geom_line(aes(x= sampling.year, y=mean.titer.obs, group = age.class, color= age.class ))+
    geom_point(aes(x= sampling.year, y=mean.titer.obs, group = age.class, color= age.class, size=N ))+
    facet_wrap(vars(age.class))+
    theme_bw()+
    ylab("log2 titer") +
    xlab('')+
    ylim(c(1,4.8)) +
    scale_x_continuous(breaks=seq(min.year.sampling,max.year.sampling,by=2))+
    theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
          axis.text.y = element_text(size=16),
          text=element_text(size=16),
          strip.background = element_blank(),
          strip.text.x = element_blank())
  return(d2)


}


# Data : Seroprevalence by age group and sampling year ----

plot_seroprevalence_age_group<- function(data, titer.threshold = 1){

  d1 =  get_data_seroprevalence(data,titer.threshold = titer.threshold)%>%
    mutate(age.class =  case_when(age <=3 ~ "1-3",
                                  age >=4 & age <= 6 ~ "4-6",
                                  age> 6 ~ "7-12" ) )%>%
    group_by(age.class, sampling.year) %>%
    mutate(NS = sum(n.samples))%>%
    mutate(Serop = sum(seropositive)) %>%
    select(sampling.year,age.class, NS,Serop)%>%
    distinct() %>%
    mutate(seroprevalence = Serop/NS) %>%
    ggplot()+
    geom_line(aes(x= sampling.year, y=seroprevalence, group = age.class, color= age.class ))+
    geom_point(aes(x= sampling.year, y=seroprevalence, group = age.class, color= age.class, size=NS))+
    theme_bw()+
    ylab("Seroprevalence") +
    xlab('')+
    ylim(c(0,1))+
    scale_x_continuous(breaks=seq(min.year.sampling,max.year.sampling))+
    theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))

  return(d1)

}

## Simulations of the seroprevalence by age group and sampling year ----

simulate_seroprevalence_age_group <- function(k, results){
  par <- results$params[k,]
  all.params = results$model$get_all_parameters(params = par,
                                                fct_model_antibody_increase = results$model$fct_model_antibody_increase,
                                                fct_model_antibody_decrease = results$model$fct_model_antibody_decrease)

  td = all.params$titer.distribution

  A = left_join(td, results$data, by = c("birth.year", "age", "titer.class", "sampling.year")) %>%
    mutate(age.class =  case_when(age <=3 ~ "1-3",
                                  age >=4 & age <= 6 ~ "4-6",
                                  age> 6 ~ "7-12" ) )%>%
    group_by(age.class, sampling.year) %>%
    mutate(N = sum(n)) %>%
    mutate(simul.titer = rmultinom(n = 1, size = N,prob = obs.proportion) ) %>%
    mutate(N.positive = sum(simul.titer[titer.class > titer.threshold])) %>% # sum on the rows that correspond to titer >titer.threshold
    mutate(seroprevalence = N.positive/N)%>%
    select(sampling.year, age.class, N, N.positive,seroprevalence) %>%
    distinct()%>%
    bind_rows()

  return(A)

}

## !!!!!!!!!!!!!!!!!! ----
## Mettre titer.threhold comme argument de la fonction simulate_seroprevalence_age_group
## !!!!!!!!!!!!!!!!!!

plot_seroprevalence_age_group_simulation <- function(results, titer.threshold = 1,burn_in = 5000, n.sim = 100){
  data = results$data

  data.sero  =  get_data_seroprevalence(data,titer.threshold = titer.threshold)%>%
    mutate(age.class =  case_when(age <=3 ~ "1-3",
                                  age >=4 & age <= 6 ~ "4-6",
                                  age> 6 ~ "7-12" ) )%>%
    group_by(age.class, sampling.year) %>%
    mutate(NS = sum(n.samples))%>%
    mutate(Serop = sum(seropositive)) %>%
    select(sampling.year,age.class, NS,Serop)%>%
    distinct() %>%
    mutate(seroprevalence = Serop/NS)

  params = results$params[-seq(1,burn_in),]

  indices = sample(x = seq(1,nrow(params)), n.sim)
  d2=indices %>%
    map(simulate_seroprevalence_age_group, results= results) %>%
    bind_rows() %>%
    group_by(age.class, sampling.year) %>%
    summarise_at(.vars = "seroprevalence",
                 .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975"))%>%
    left_join(data.sero , by = c( "age.class", "sampling.year")) %>%
    ggplot()+
    geom_line(aes(x= sampling.year, y= mean, group = age.class, color = age.class), linetype = 'dashed')+
    geom_ribbon(aes(x = sampling.year, ymin=quantile025, ymax=quantile975,group = age.class,fill=age.class), alpha=0.2) +
    geom_line(aes(x= sampling.year, y=seroprevalence, group = age.class, color= age.class ))+
    geom_point(aes(x= sampling.year, y=seroprevalence, group = age.class, color= age.class, size=NS))+
    facet_wrap(vars(age.class))+
    theme_bw()+
    ylab("Seroprevalence") +
    xlab('')+
    ylim(c(0,1)) +
    scale_x_continuous(breaks=seq(min.year.sampling,max.year.sampling, by=2))+
    theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
          axis.text.y = element_text(size=16),
          text=element_text(size=16),
          strip.background = element_blank(),
          strip.text.x = element_blank())

  return(d2)

}

