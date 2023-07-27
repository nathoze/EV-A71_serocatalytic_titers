# Definition of the function specific for the models
# Names of the functions : get_all_parameters_model1,  update_all_parameters_model1 etc
# to replace the generic functions  get_all_parameters, update_all_parameters

# list the functions to be changed here


# Only one parameter for the FOI in the constant model
update_all_parameters_model_constant <- function(old.all.parameters, new.param, updated_index, fct_model_antibody_increase, fct_model_antibody_decrease){

  new.all.parameters = old.all.parameters
  new.all.parameters$params[updated_index] = new.param

  if(updated_index == 1){ # We change the foi and therefore only the corresponding matrix
    p = infection_probability(new.param)
    for(j in 1:N.FOI){
      new.all.parameters$transition.matrices[[j]] = (1-p)*new.all.parameters$decay.matrix+p*new.all.parameters$increase.matrix
    }
  #  new.all.parameters$transition.matrices[[updated_index]] = (1-p)*new.all.parameters$decay.matrix+p*new.all.parameters$increase.matrix
  }

  if(updated_index >= 2){ #  Increase  or decay parameter
    new.all.parameters = get_all_parameters_model_constant(new.all.parameters$params, fct_model_antibody_increase, fct_model_antibody_decrease)
  }

  new.all.parameters$titer.distribution = titer_distribution(new.all.parameters$transition.matrices)

  return(new.all.parameters)

}

get_all_parameters_model_constant <- function(params, fct_model_antibody_increase, fct_model_antibody_decrease){
  p = infection_probability(params[1])
  increase.matrix=fct_model_antibody_increase(N.titers  = N.titers, sigma.P = params[2])
  decay.matrix=fct_model_antibody_decrease(N.titers = N.titers, omega = params[3])

  transition.matrices = NULL
  for(i in 1:N.FOI){
     transition.matrices[[i]] = (1-p)*decay.matrix+p*increase.matrix
  }

  titer.distribution = titer_distribution(transition.matrices)

  all.parameters = list(increase.matrix = increase.matrix,
                        decay.matrix = decay.matrix,
                        transition.matrices = transition.matrices,
                        titer.distribution = titer.distribution,
                        params = params)


  return(all.parameters)

}

is_invalid_model_constant <- function(k, value) { # Function that checks if parameter value is invalid
  if (value <10e-9) { return(TRUE) } # All the parameters must be > 0
  if (k == 1 & value >2) { return(TRUE) } # the foi
  if (k == 2 & value >8) { return(TRUE) }# sigmaP
  if (k == 3 & value >8) { return(TRUE) }# Omega
  FALSE
}


# Model where the FOI is constant during five years

update_all_parameters_model_five_years <- function(old.all.parameters, new.param, updated_index, fct_model_antibody_increase, fct_model_antibody_decrease){

  new.all.parameters = old.all.parameters
  new.all.parameters$params[updated_index] = new.param

  if(updated_index <= round(N.FOI/5)){ # We change the foi and therefore only the corresponding matrix
    p = infection_probability(new.param)
    #J =   floor((i-1)/5) + 1
    for(k in seq( (updated_index-1)*5+1,(updated_index)*5) ){
      new.all.parameters$transition.matrices[[k]] = (1-p)*new.all.parameters$decay.matrix+p*new.all.parameters$increase.matrix
    }
 #   new.all.parameters$transition.matrices[[updated_index]] = (1-p)*new.all.parameters$decay.matrix+p*new.all.parameters$increase.matrix
  }

  if(updated_index > round(N.FOI/5)){ #  Increase  or decay parameter
    new.all.parameters = get_all_parameters_model_five_years(new.all.parameters$params, fct_model_antibody_increase, fct_model_antibody_decrease)
  }

  new.all.parameters$titer.distribution = titer_distribution(new.all.parameters$transition.matrices)

  return(new.all.parameters)

}
get_all_parameters_model_five_years <- function(params, fct_model_antibody_increase, fct_model_antibody_decrease){

  increase.matrix=fct_model_antibody_increase(N.titers  = N.titers, sigma.P = head(tail(params, n=2), n=1))
  decay.matrix=fct_model_antibody_decrease(N.titers = N.titers, omega = tail(params,1))

  transition.matrices = NULL
  for(i in 1:N.FOI){
    J =   floor((i-1)/5) + 1
    p = infection_probability(params[J])
    transition.matrices[[i]] = (1-p)*decay.matrix+p*increase.matrix
  }

  titer.distribution = titer_distribution(transition.matrices)

  all.parameters = list(increase.matrix = increase.matrix,
                        decay.matrix = decay.matrix,
                        transition.matrices = transition.matrices,
                        titer.distribution = titer.distribution,
                        params = params)


  return(all.parameters)

}
is_invalid_model_five_years <- function(k, value) { # Function that checks if parameter value is invalid
  if (value <10e-9) { return(TRUE) } # All the parameters must be > 0
  if (k<=6  & value >2) { return(TRUE) } # the foi
  if (k == 7 & value >8) { return(TRUE) }# sigmaP
  if (k == 8 & value >8) { return(TRUE) }# Omega
  FALSE
}


get_all_parameters_model_independent <- function(params, fct_model_antibody_increase, fct_model_antibody_decrease){

  increase.matrix=get_increase_matrix(N.titers  = N.titers, sigma.P = params[N.FOI+1])
  decay.matrix=get_decay_matrix(N.titers = N.titers, omega = params[N.FOI+2])
  transition.matrices = NULL

  for(i in 1:N.FOI){
    p = infection_probability(params[i])
    transition.matrices[[i]] = (1-p)*decay.matrix+p*increase.matrix
  }

  titer.distribution = titer_distribution(transition.matrices)

  all.parameters = list(increase.matrix = increase.matrix,
                        decay.matrix = decay.matrix,
                        transition.matrices = transition.matrices,
                        titer.distribution = titer.distribution,
                        params = params)


  return(all.parameters)

}

# When the parameters are updated in the MCMC, only update some of the transformed parameters
update_all_parameters_model_independent <- function(old.all.parameters, new.param, updated_index, fct_model_antibody_increase, fct_model_antibody_decrease){

  new.all.parameters = old.all.parameters
  new.all.parameters$params[updated_index] = new.param

  if(updated_index <= N.FOI){ # We change the foi and therefore only the corresponding matrix
    p = infection_probability(new.param)
    new.all.parameters$transition.matrices[[updated_index]] = (1-p)*new.all.parameters$decay.matrix+p*new.all.parameters$increase.matrix

    # I "optimized" the update of parameters but it didn't change the runtime
    # t.d = titer_distribution_optimized(new.all.parameters, updated_index)
    # new.all.parameters$titer.distribution = new.all.parameters$titer.distribution %>%
    #   left_join(t.d, by = c("birth.year", "age", "titer.class", "sampling.year")) %>%
    #   mutate(obs.proportion = case_when(is.na(obs.proportion.y) ~ obs.proportion.x, TRUE ~ obs.proportion.y)) %>%
    #   select(-c(obs.proportion.x, obs.proportion.y))
  }

  if(updated_index >= N.FOI+1){ #  Increase  or decay parameter
    new.all.parameters = get_all_parameters_model_independent(new.all.parameters$params, fct_model_antibody_increase, fct_model_antibody_decrease)
  }

  new.all.parameters$titer.distribution = titer_distribution(new.all.parameters$transition.matrices)

  return(new.all.parameters)

}
is_invalid_model_independent <- function(k, value) { # Function that checks if parameter value is invalid
  if (value <10e-9) { return(TRUE) } # All the parameters must be > 0
  if (k<=N.FOI & value >2) { return(TRUE) } # the foi
  if (k >= N.FOI+1 & value >8) { return(TRUE) }# sigmaP, Omega
  #if (k == N.FOI+2 & value >8) { return(TRUE) }# Omega
  FALSE
}

