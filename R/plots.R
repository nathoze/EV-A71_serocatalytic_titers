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
    facet_wrap(  vars(sampling.year))+
    scale_x_continuous(breaks=seq(1,12))+
    ylim(c(0, 5))+
    ylab('log2 Titer')+
    theme_bw()+
    facet_wrap(vars(sampling.year))

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
    theme_bw()+
 #   #theme(axis.text.x = element_text(size=14),
      #                axis.text.y = element_text(  size=14),
     #                 text=element_text(size=14))+
    facet_wrap(vars(sampling.year))

  return(g)

}
