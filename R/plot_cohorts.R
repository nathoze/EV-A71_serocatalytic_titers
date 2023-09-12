library(gridExtra)
# plot cohorts -----
results =  readRDS(file='results/Model_independent_no_protection.rds')
data = results$data
burn_in = 5000
n.sim=100
params = results$params[-seq(1,burn_in),]
indices = sample(x = seq(1,nrow(params)), n.sim)
## titer vs sampling year ----

g1 = data%>%
  group_by(age,birth.year) %>%
  mutate(N = sum(n))    %>%
  mutate(birth.year = factor(birth.year, levels = as.character(birth.years))) %>%
  mutate(mean.titer.obs = sum(titer.class*n)/N) %>%
  group_by(sampling.year) %>%
  mutate(N2 = sum(n))    %>%
  mutate(mean.titer.by.year  = sum(titer.class*n)/N2 ) %>%
  group_by(age,birth.year) %>%
  ggplot()+
  geom_line(aes(x=sampling.year, y = mean.titer.obs, group = birth.year, color= birth.year)) +
  geom_point(aes(x=sampling.year, y = mean.titer.obs, group = birth.year, color= birth.year)) +
  geom_line(aes(x=sampling.year, y = mean.titer.by.year), color= "black",linewidth=1.2) +
  geom_point(aes(x=sampling.year, y = mean.titer.by.year),   color= "black",size=1.6) +
  scale_color_discrete(name="Cohort")+
  ylab('log2 Titer')+
  ylim(c(0, 5))+
  theme_bw() +
  xlab('Sampling year')+
  scale_x_continuous(breaks=seq(1995,2012))+
  theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
        axis.text.y = element_text(size=16),
        text=element_text(size=16))+
  ylim(c(0,NA))
print(g1)

dev.copy(pdf,"results/cohort_titer_year.pdf", width =9, height = 6)
dev.off()

## titer vs age ----

g2= data%>%
  group_by(age,birth.year) %>%
  mutate(N = sum(n))    %>%
  mutate(birth.year = factor(birth.year, levels = as.character(birth.years))) %>%
  mutate(mean.titer.obs = sum(titer.class*n)/N) %>%
  # ungroup()%>%
  group_by(age)%>%
  mutate(N2 = sum(n))    %>%
  mutate(mean.titer.age = sum(titer.class*n)/N2) %>%
  #ungroup()%>%
  group_by(age,birth.year) %>%
  ggplot()+
  geom_line(aes(x=age, y = mean.titer.obs, group = birth.year, color= birth.year)) +
  geom_point(aes(x=age, y = mean.titer.obs, group = birth.year, color= birth.year)) +
  geom_line(aes(x=age, y = mean.titer.age), color= "black",linewidth=1.2) +
  geom_point(aes(x=age, y = mean.titer.age),   color= "black",size=1.6) +
  scale_color_discrete(name="Cohort")+
  ylab('log2 Titer')+
  ylim(c(0, 5))+
  theme_bw() +
  xlab('Age')+
  scale_x_continuous(breaks=seq(1,12))+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        text=element_text(size=16))+
  ylim(c(0,NA))
print(g2)

dev.copy(pdf,"results/cohort_titer_age.pdf", width =9, height = 6)
dev.off()

titer.threshold = 2
if(titer.threshold ==1){
  data.seroprevalence = data_1
}
if(titer.threshold ==2){
  data.seroprevalence = data_2
}
if(titer.threshold ==3){
  data.seroprevalence = data_3
}

## Seroprevalence vs Sampling year ----

g3 = data.seroprevalence%>%
  group_by(age,birth.year) %>%
  # mutate(N = sum(n))    %>%
  mutate(birth.year = factor(birth.year, levels = as.character(birth.years))) %>%
  ungroup()%>%
  group_by(sampling.year)%>%
  mutate(mean.seroprevalence = sum(seropositive)/sum(n.samples)) %>%
  ungroup()%>%
  group_by(age,birth.year) %>%
  #mutate(mean.titer.obs = sum(titer.class*n)/N) %>%
  ggplot()+
  geom_line(aes(x=sampling.year, y = seroprevalence, group = birth.year, color= birth.year)) +
  geom_point(aes(x=sampling.year, y = seroprevalence, group = birth.year, color= birth.year)) +
  geom_line(aes(x=sampling.year, y = mean.seroprevalence ), color= "black",linewidth=1.2) +
  geom_point(aes(x=sampling.year, y = mean.seroprevalence),   color= "black",size=1.6) +
  scale_color_discrete(name="Cohort")+
  ylab('Seroprevalence')+
  ylim(c(0,1))+
  theme_bw() +
  xlab('Sampling year')+
  scale_x_continuous(breaks=seq(1995,2012))+
  theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
        axis.text.y = element_text(size=16),
        text=element_text(size=16))
print(g3)

dev.copy(pdf,"results/cohort_seroprevalence_year.pdf", width =9, height = 6)
dev.off()

## Seroprevalence vs age ----

g4 = data.seroprevalence%>%
  group_by(age,birth.year) %>%
  mutate(birth.year = factor(birth.year, levels = as.character(birth.years))) %>%
  group_by(age) %>%
  mutate(mean.seroprevalence = sum(seropositive)/sum(n.samples)) %>%
  group_by(age,birth.year) %>%
  ggplot()+
  geom_line(aes(x = age, y = seroprevalence, group = birth.year, color= birth.year)) +
  geom_point(aes(x = age, y = seroprevalence, group = birth.year, color= birth.year)) +
  geom_line(aes(x = age, y = mean.seroprevalence, group =birth.year  ), color= "black",linewidth=1.2) +
  geom_point(aes(x = age, y = mean.seroprevalence, group =birth.year),   color= "black",size=1.6) +
  scale_color_discrete(name="Cohort")+
  ylab('Seroprevalence')+
  ylim(c(0,1))+
  theme_bw() +
  xlab('Age')+
  scale_x_continuous(breaks=seq(12))+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        text=element_text(size=16))
print(g4)

dev.copy(pdf,"results/cohort_seroprevalence_age.pdf", width =9, height = 6)
dev.off()

###################
## add the fits ---
##################

## titer vs sampling year ----

data.plot = data%>%
  group_by(age,birth.year) %>%
  mutate(N = sum(n))    %>%
  mutate(birth.year = factor(birth.year, levels = as.character(birth.years))) %>%
  mutate(mean.titer.obs = sum(titer.class*n)/N) %>%
  group_by(sampling.year) %>%
  mutate(N2 = sum(n))    %>%
  mutate(mean.titer.by.year  = sum(titer.class*n)/N2 ) %>%
  group_by(age,birth.year)




simulate_mean_posterior_titer_sampling_year <-function(k, results){

  par <- results$params[k,]
  all.params = results$model$get_all_parameters(params = par,
                                                fct_model_antibody_increase = results$model$fct_model_antibody_increase,
                                                fct_model_antibody_decrease = results$model$fct_model_antibody_decrease)

  td = all.params$titer.distribution

  A = left_join(td, results$data, by = c("birth.year", "age", "titer.class", "sampling.year")) %>%
    group_by(sampling.year) %>%
    mutate(N = sum(n)) %>%
    mutate(simul.titer = rmultinom(n = 1, size = N,prob = obs.proportion)*titer.class/N ) %>%
    summarise(mean.titer.sim=colSums(simul.titer)) #%>% mutate(index = k)

  return(A)

}


#
# g = data.plot%>%
#   ggplot()+
#   geom_line(aes(x=sampling.year, y = mean.titer.obs, group = birth.year, color= birth.year)) +
#   geom_point(aes(x=sampling.year, y = mean.titer.obs, group = birth.year, color= birth.year)) +
#   geom_line(aes(x=sampling.year, y = mean.titer.by.year), color= "black",linewidth=1.2) +
#   geom_point(aes(x=sampling.year, y = mean.titer.by.year),   color= "black",size=1.6) +
#   scale_color_discrete(name="Cohort")+
#   ylab('log2 Titer')+
#   ylim(c(0, 5))+
#   theme_bw() +
#   xlab('Sampling year')+
#   scale_x_continuous(breaks=seq(1995,2012))+
#   theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
#         axis.text.y = element_text(size=16),
#         text=element_text(size=16))+
#   ylim(c(0,NA))
# print(g)




g5 = indices %>%
  map(simulate_mean_posterior_titer_sampling_year, results= results) %>%
  bind_rows() %>%
  group_by(sampling.year)%>%
  summarise_at(.vars = "mean.titer.sim",
               .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
  left_join(data.plot , by = c( "sampling.year")) %>%
  ggplot()+
  geom_line(aes(x= sampling.year, y= mean), linetype = 'dashed', color= "red")+
  geom_ribbon(aes(x = sampling.year, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= "red") +
  geom_line(aes(x=sampling.year, y = mean.titer.obs, group = birth.year, color= birth.year)) +
  geom_point(aes(x=sampling.year, y = mean.titer.obs, group = birth.year, color= birth.year)) +
  geom_line(aes(x=sampling.year, y = mean.titer.by.year), color= "black",linewidth=1.2) +
  geom_point(aes(x=sampling.year, y = mean.titer.by.year),   color= "black",size=1.6) +
  scale_color_discrete(name="Cohort")+
  ylab('log2 Titer')+
  ylim(c(0, 5))+
  theme_bw() +
  xlab('Sampling year')+
  scale_x_continuous(breaks=seq(1995,2012))+
  theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
        axis.text.y = element_text(size=16),
        text=element_text(size=16))+
  ylim(c(0,NA))
print(g5)


dev.copy(pdf,"results/cohort_titer_year_sim.pdf", width =9, height = 6)
dev.off()




## titer vs age ----

data.plot = data%>%
  group_by(age,birth.year) %>%
  mutate(N = sum(n))    %>%
  mutate(birth.year = factor(birth.year, levels = as.character(birth.years))) %>%
  mutate(mean.titer.obs = sum(titer.class*n)/N) %>%
  group_by(age)%>%
  mutate(N2 = sum(n))    %>%
  mutate(mean.titer.age = sum(titer.class*n)/N2) %>%
  group_by(age,birth.year)


simulate_mean_posterior_titer_all_ages <-function(k, results){

  par <- results$params[k,]
  all.params = results$model$get_all_parameters(params = par,
                                                fct_model_antibody_increase = results$model$fct_model_antibody_increase,
                                                fct_model_antibody_decrease = results$model$fct_model_antibody_decrease)

  td = all.params$titer.distribution

  A = left_join(td, results$data, by = c("birth.year", "age", "titer.class", "sampling.year")) %>%
    group_by(age) %>%
    mutate(N = sum(n)) %>%
    mutate(simul.titer = rmultinom(n = 1, size = N,prob = obs.proportion)*titer.class/N ) %>%
    summarise(mean.titer.sim=colSums(simul.titer)) #%>% mutate(index = k)

  return(A)

}


g6 = indices %>%
  map(simulate_mean_posterior_titer_all_ages, results= results) %>%
  bind_rows() %>%
  group_by(age)%>%
  summarise_at(.vars = "mean.titer.sim",
               .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
  left_join(data.plot , by = c( "age")) %>%
  ggplot()+
  geom_line(aes(x= age, y= mean), linetype = 'dashed', color= "red")+
  geom_ribbon(aes(x = age, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= "red") +
  geom_line(aes(x=age, y = mean.titer.obs, group = birth.year, color= birth.year)) +
  geom_point(aes(x=age, y = mean.titer.obs, group = birth.year, color= birth.year)) +
  geom_line(aes(x=age, y = mean.titer.age), color= "black",linewidth=1.2) +
  geom_point(aes(x=age, y = mean.titer.age),   color= "black",size=1.6) +
  scale_color_discrete(name="Cohort")+
  ylab('log2 Titer')+
  ylim(c(0, 5))+
  theme_bw() +
  xlab('Age')+
  scale_x_continuous(breaks=seq(1,12))+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        text=element_text(size=16))+
  ylim(c(0,NA))

print(g6)


dev.copy(pdf,"results/cohort_titer_age_sim.pdf", width =9, height = 6)
dev.off()

## Seroprevalence vs Sampling year ----


data.plot = data.seroprevalence%>%
  group_by(age,birth.year) %>%
  # mutate(N = sum(n))    %>%
  mutate(birth.year = factor(birth.year, levels = as.character(birth.years))) %>%
  ungroup()%>%
  group_by(sampling.year)%>%
  mutate(mean.seroprevalence = sum(seropositive)/sum(n.samples)) %>%
  ungroup()%>%
  group_by(age,birth.year)




simulate_seroprevalence_sampling_year <- function(k, results){
  par <- results$params[k,]
  all.params = results$model$get_all_parameters(params = par,
                                                fct_model_antibody_increase = results$model$fct_model_antibody_increase,
                                                fct_model_antibody_decrease = results$model$fct_model_antibody_decrease)

  td = all.params$titer.distribution

  A = left_join(td, results$data, by = c("birth.year", "age", "titer.class", "sampling.year")) %>%
    group_by(sampling.year) %>%
    mutate(N = sum(n)) %>%
    select(age,obs.proportion,titer.class,N) %>%
    distinct() %>%
    mutate(simul.titer = rmultinom(n = 1, size = N,prob = obs.proportion) ) %>%
    mutate(N.positive = sum(simul.titer[titer.class > titer.threshold])) %>% # sum on the rows that correspond to titer >titer.threshold
    mutate(seroprevalence = N.positive/N)%>%
    select(age,seroprevalence) %>%
    distinct()%>%
    bind_rows()

  return(A)

}

g7 = indices %>%
  map(simulate_seroprevalence_sampling_year, results= results) %>%
  bind_rows() %>%
  summarise_at(.vars = "seroprevalence",
               .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
  left_join(data.plot , by = c( "sampling.year")) %>%
  ggplot()+
  geom_line(aes(x= sampling.year, y= mean), linetype = 'dashed', color= "red")+
  geom_ribbon(aes(x = sampling.year, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= "red") +
  geom_line(aes(x=sampling.year, y = seroprevalence, group = birth.year, color= birth.year)) +
  geom_point(aes(x=sampling.year, y = seroprevalence, group = birth.year, color= birth.year)) +
  geom_line(aes(x=sampling.year, y = mean.seroprevalence ), color= "black",linewidth=1.2) +
  geom_point(aes(x=sampling.year, y = mean.seroprevalence),   color= "black",size=1.6) +
  scale_color_discrete(name="Cohort")+
  ylab('Seroprevalence')+
  ylim(c(0,1))+
  theme_bw() +
  xlab('Sampling year')+
  scale_x_continuous(breaks=seq(1995,2012))+
  theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
        axis.text.y = element_text(size=16),
        text=element_text(size=16))
print(g7)

dev.copy(pdf,"results/cohort_seroprevalence_year_sim.pdf", width =9, height = 6)
dev.off()


## Seroprevalence vs age ----



data.plot  = data.seroprevalence%>%
  group_by(age,birth.year) %>%
  mutate(birth.year = factor(birth.year, levels = as.character(birth.years))) %>%
  group_by(age) %>%
  mutate(mean.seroprevalence = sum(seropositive)/sum(n.samples)) %>%
  group_by(age,birth.year)



# titer.threshold en argument ??
simulate_seroprevalence_all_ages <- function(k, results){
  par <- results$params[k,]
  all.params = results$model$get_all_parameters(params = par,
                                                fct_model_antibody_increase = results$model$fct_model_antibody_increase,
                                                fct_model_antibody_decrease = results$model$fct_model_antibody_decrease)

  td = all.params$titer.distribution

  A = left_join(td, results$data, by = c("birth.year", "age", "titer.class", "sampling.year")) %>%
    group_by(age) %>%
    mutate(N = sum(n)) %>%
    select(age,obs.proportion,titer.class,N) %>%
    distinct() %>%
    mutate(simul.titer = rmultinom(n = 1, size = N,prob = obs.proportion) ) %>%
    mutate(N.positive = sum(simul.titer[titer.class > titer.threshold])) %>% # sum on the rows that correspond to titer >titer.threshold
    mutate(seroprevalence = N.positive/N)%>%
    select(age,seroprevalence) %>%
    distinct()%>%
    bind_rows()

  return(A)

}



g8 = indices %>%
  map(simulate_seroprevalence_all_ages, results= results) %>%
  bind_rows() %>%
  summarise_at(.vars = "seroprevalence",
               .funs = c(mean="mean",quantile025 = "quantile025", quantile975="quantile975")) %>%
  left_join(data.plot , by = c( "age")) %>%
  ggplot()+
  geom_line(aes(x= age, y= mean), linetype = 'dashed', color= "red")+
  geom_ribbon(aes(x = age, ymin=quantile025, ymax=quantile975), alpha=0.2, fill= "red") +
  geom_line(aes(x = age, y = seroprevalence, group = birth.year, color= birth.year)) +
  geom_point(aes(x = age, y = seroprevalence, group = birth.year, color= birth.year)) +
  geom_line(aes(x = age, y = mean.seroprevalence, group =birth.year  ), color= "black",linewidth=1.2) +
  geom_point(aes(x = age, y = mean.seroprevalence, group =birth.year),   color= "black",size=1.6) +
  scale_color_discrete(name="Cohort")+
  ylab('Seroprevalence')+
  ylim(c(0,1))+
  theme_bw() +
  xlab('Age')+
  scale_x_continuous(breaks=seq(12))+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        text=element_text(size=16))
print(g8)


dev.copy(pdf,"results/cohort_seroprevalence_age_sim.pdf", width =9, height = 6)
dev.off()











#
# g = data.seroprevalence%>%
#   group_by(age,birth.year) %>%
#   # mutate(N = sum(n))    %>%
#   mutate(birth.year = factor(birth.year, levels = as.character(birth.years))) %>%
#   #mutate(mean.titer.obs = sum(titer.class*n)/N) %>%
#   ggplot()+
#   geom_line(aes(x=sampling.year, y = seroprevalence, group = birth.year, color= birth.year)) +
#   geom_point(aes(x=sampling.year, y = seroprevalence, group = birth.year, color= birth.year)) +
#   facet_wrap(vars(birth.year))+
#   scale_color_discrete(name="Cohort")+
#   ylab('Seroprevalence')+
#   ylim(c(0,1))+
#   theme_bw() +
#   xlab('Sampling year')+
#   scale_x_continuous(breaks=seq(1995,2012,by=4))+
#   theme(axis.text.x = element_text(size=16,angle = 45, vjust=0.5),
#         axis.text.y = element_text(size=16),
#         text=element_text(size=16))
# print(g)
#

grid.arrange(g1,g2,g3,g4,
             ncol =2, nrow = 2)

dev.copy(pdf,"results/cohort_data.pdf", width=18, height =12)
dev.off()

grid.arrange(g5,g6,g7,g8,
             ncol =2, nrow = 2)

dev.copy(pdf,"results/cohort_sim.pdf", width=18, height =12)
dev.off()

