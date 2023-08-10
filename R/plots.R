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

plot_fit <- function(results, burn_in, n.sim=50){

  data = results$data%>%
    group_by(age, sampling.year) %>%
    mutate(N = sum(n)) %>%
    mutate(mean.titer.obs = sum(titer.class*n)/N) %>%
    select(sampling.year,age , mean.titer.obs )

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

  return(g)

}

plot_data <- function(data){
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
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          text=element_text(size=16))+
    ylim(c(0,NA))

  return(g)

}
