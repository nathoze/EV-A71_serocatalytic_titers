# plot the fits of the model to the data
# simulate data from the posterior distribution

plot_fit <- function(results, burn_in){

  data = results$data
  params = results$params[-seq(1,burn_in),]
  n.sim = 100

   # --> summarise data : get total number by age group and sampling year

  #indices = sample(x = seq(1,nrow(params)), n.sim)
  par <- colMeans(params)

  all.params = results$model$get_all_parameters(params = par,
                                                fct_model_antibody_increase = results$model$fct_model_antibody_increase,
                                                fct_model_antibody_decrease = results$model$fct_model_antibody_decrease)

  td= all.params$titer.distribution

  A = left_join(td, data, by = c("birth.year", "age", "titer.class", "sampling.year")) %>%
    group_by(age, sampling.year) %>%
    mutate(N = sum(n)) %>%
    mutate(simul.titer = rmultinom(n = n.sim, size = N,prob = obs.proportion)*titer.class/N ) %>%
    summarise(mean.titer.sim=colSums(simul.titer))

  data2=data %>%
    group_by(age, sampling.year) %>%
    mutate(N = sum(n)) %>%
    mutate(mean.titer.obs = sum(titer.class*n)/N)   %>%
    select(sampling.year,age , mean.titer.obs )


  g= A %>%
    left_join(data2, by = c( "age", "sampling.year")) %>%
    ggplot()+
    stat_summary(aes(x=age, y = mean.titer.sim),
                 fun.data=mean_sdl, fun.args = list(mult=1),
                 geom="pointrange", color="red")+
    geom_line(aes(x=age, y = mean.titer.obs)) +
    geom_point(aes(x=age, y = mean.titer.obs)) +
    facet_wrap(  vars(sampling.year))+
    scale_x_continuous(breaks=seq(1,12))+
    theme_bw()+
    facet_wrap(  vars(sampling.year))

  return(g)



}
