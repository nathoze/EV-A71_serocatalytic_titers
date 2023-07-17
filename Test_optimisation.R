

get_increase_matrix <- function(N.titers = 10, sigma.P = 1){ 
  
  Titers <- seq(1,N.titers)
  Matrix_transition = matrix(data=NA, nrow = N.titers, ncol = N.titers)
  for(titer in Titers){
    if(titer < N.titers){ 
      infection.pmf = truncated_poisson_pmf(sigma.P,truncation.point = N.titers-titer,Titers)
      infection.pmf = c(rep(0, titer), head(infection.pmf,N.titers-titer)) #[-1]
    } 
    Matrix_transition[titer,] = infection.pmf
    
  }
  return(Matrix_transition)
}
 
get_transition_matrix <- function(foi, N.titers = 10, sigma.P = 1, omega=1 ){ 
  
  Titers <- seq(1,N.titers)
  
  proba.infection = infection_probability(foi)
  Matrix_transition = matrix(data=NA, nrow = N.titers, ncol = N.titers)
  for(titer in Titers){
    if(titer < N.titers){ 
      infection.pmf = truncated_poisson_pmf(sigma.P,truncation.point = N.titers-titer,Titers)
      infection.pmf = head(infection.pmf,N.titers-titer)#[-1]
    } 
    
    if(titer == min(Titers) ){ 
      decay = 1
    } 
    
    if(titer > min(Titers) ){ 
      # decay =rep(0,N.titers)
      decay = decay_poisson_pmf(omega,truncation.point = titer-1,0:N.titers)
      decay = tail(decay,titer)
    }
    
    A=c(decay*(1-proba.infection),proba.infection*infection.pmf)
    if(titer == N.titers){ 
      A = decay#*(1-proba.infection)
    }
    
    Matrix_transition[titer,] = A
    
  }
  return(Matrix_transition)
}

get_decay_matrix <- function(N.titers = 10,omega = 1){ 
  
  Titers <- seq(1,N.titers)
  Matrix_transition = matrix(data=NA, nrow = N.titers, ncol = N.titers)
  for(titer in Titers){
    
    if(titer == min(Titers) ){ 
      decay = 1
    } 
    if(titer > min(Titers) ){ 
      # decay =rep(0,N.titers)
      decay = decay_poisson_pmf(omega,truncation.point = titer-1,0:N.titers)
      decay = tail(decay,titer)
    }
    
    A=c(decay, rep(0, N.titers-titer))
    
    if(titer == N.titers){ 
      A = decay
    }
    Matrix_transition[titer,] = A
  }
  
  return(Matrix_transition)
}

infection_probability_2 <- function(foi,N.titers = 10){
  return(c(rep(1-exp(-foi),N.titers-1), 0))
}


foi= 0.6
p = infection_probability_2(foi)

omega=1
sigma.P=2
A1= (1-p)*get_decay_matrix(10,omega)+p*get_increase_matrix(N.titers  = 10, sigma.P = sigma.P)
A2 = get_transition_matrix(foi,10,sigma.P = sigma.P,omega = omega)

# test that both matrices are equal
norm(A1-A2)



foi.list = data.frame(FOI.years = birth.years, FOI = lambda )


transition.matrices = NULL
i=0
for(foi in foi.list$FOI){  
  i=i+1
  p = infection_probability_2(foi)
  transition.matrices[[i]] = (1-p)*get_decay_matrix(10,omega)+p*get_increase_matrix(N.titers  = 10, sigma.P = sigma.P)
}



start_time <- Sys.time()
transition.matrices = NULL
i=0
for(foi in foi.list$FOI){  
  i=i+1
  p = infection_probability_2(foi)
  transition.matrices[[i]] = (1-p)*get_decay_matrix(10,omega)+p*get_increase_matrix(N.titers  = 10, sigma.P = sigma.P)
}
end_time <- Sys.time()
print(end_time - start_time)


start_time <- Sys.time()
decay.matrix=get_decay_matrix(10,omega)
increase.matrix=get_increase_matrix(N.titers  = 10, sigma.P = sigma.P)
transition.matrices = NULL
i=0
for(foi in foi.list$FOI){  
  i=i+1
  p = infection_probability_2(foi)
  transition.matrices[[i]] = (1-p)*decay.matrix+p*increase.matrix
}
end_time <- Sys.time()
print(end_time - start_time)

# Separate the computation of the titers and the computation of the likelihood

# All.parameters : parameters + transformed parameters (Matrices + titer distribution)

compute_loglik_2 <- function(All.titers) {
  
  prob.titers =  left_join(All.titers, data.EV71.Malaysia,
                           by = c("birth.year", "age", "titer.class", "sampling.year"))  %>% 
    group_by(age, sampling.year) %>% 
    summarise(ll  = dmultinom(x = n,prob = obs.proportion, log=TRUE), .groups='drop')  %>%
    summarise(s= sum(ll))
  
  return(prob.titers$s)
  
}

 
get_all_probability_titers_2 <- function(all.parameters){

  df = NULL
  j=0
  for(birth.year in seq(min.year.sampling-age.max,max.year.sampling )){ 
    titers = Titers.0 
    
    ## Attention aux années prises pour le dernier échantillonage
    years <- head(seq(birth.year,min(birth.year+age.max,max.year.sampling)),age.max) 
 
    n = length(years)
    Titers.year = matrix(data = NA, nrow = N.titers, ncol = n)
    i=0
    for(foi.index in years){
      i=i+1
      titers = t(all.parameters$transition.matrices[[i+j]])%*% titers
      Titers.year[,i] = titers
    }
    j=j+1
    
    obs = get_observed_titers(Titers.year)
    
    df=rbind(df, data.frame( birth.year = birth.year, 
                             age = rep(seq(1,n), each =titer_observable_max),
                             obs.proportion = as.numeric(obs), titer.class= rep(seq(1,titer_observable_max), n)))
  }
  
  
  df2 = df %>% 
    mutate(sampling.year = age+birth.year) %>%
    #filter(age <= age.max) %>% 
    filter(sampling.year >= min.year.sampling) %>% 
    filter(sampling.year <= max.year.sampling)
  
  return(df2)
 
}


get_all_parameters <- function(params){
  
  increase.matrix=get_increase_matrix(N.titers  = N.titers, sigma.P = params[N.FOI+1])
  decay.matrix=get_decay_matrix(N.titers = N.titers, omega = params[N.FOI+2])
  transition.matrices = NULL
  
  for(i in 1:N.FOI){  
    p = infection_probability_2(params[i])
    transition.matrices[[i]] = (1-p)*decay.matrix+p*increase.matrix
  }

  all.parameters = list(increase.matrix = increase.matrix,  
                        decay.matrix = decay.matrix, 
                        transition.matrices = transition.matrices,
                        params = params)
  
  all.parameters$titer.distribution = get_all_probability_titers_2(all.parameters)
  
  return(all.parameters)
  
}

# When the parameters are updated in the MCMC, only update some of the transformed parameters
update_all_parameters <- function(old.all.parameters, new.param, updated_index){
  
  new.all.parameters = old.all.parameters
  new.all.parameters$params[updated_index] = new.param
  if(updated_index <= N.FOI){ # We change the foi and therefore only the corresponding matrix
    p = infection_probability_2(new.param)
    new.all.parameters$transition.matrices[[updated_index]] = (1-p)*new.all.parameters$decay.matrix+p*new.all.parameters$increase.matrix
  }  
  
  if(updated_index >= N.FOI+1){ #  Increase  or decay parameter
    new.all.parameters = get_all_parameters(new.all.parameters$params)
  }
  
  # TO DO : compute the titer distribution on specific indices
  new.all.parameters$titer.distribution = get_all_probability_titers_2(new.all.parameters)
  
  return(new.all.parameters)
  
}


compute_loglik(params0)
 

G=get_all_parameters(params0) 
All.titers =G$titer.distribution  
compute_loglik_2(G$titer.distribution)

