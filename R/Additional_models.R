# Definition of the function specific for the models
# Names of the functions : get_all_parameters_model1,  update_all_parameters_model1 etc
# to replace the generic functions  get_all_parameters, update_all_parameters

# list the functions to be changed here


# Only one parameter for the FOI in the constant model
update_all_parameters_model_constant <- function(old.all.parameters, new.param, updated_index){

  new.all.parameters = old.all.parameters
  new.all.parameters$params[updated_index] = new.param

  if(updated_index == 1){ # We change the foi and therefore only the corresponding matrix
    p = infection_probability(new.param)
    new.all.parameters$transition.matrices[[updated_index]] = (1-p)*new.all.parameters$decay.matrix+p*new.all.parameters$increase.matrix
  }

  if(updated_index >= 2){ #  Increase  or decay parameter
    new.all.parameters = get_all_parameters_model_constant(new.all.parameters$params)
  }

  new.all.parameters$titer.distribution = titer_distribution(new.all.parameters$transition.matrices)

  return(new.all.parameters)

}
get_all_parameters_model_constant <- function(params){

  p = infection_probability(params[1])
  increase.matrix=get_increase_matrix(N.titers  = N.titers, sigma.P = params[2])
  decay.matrix=get_decay_matrix(N.titers = N.titers, omega = params[3])

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
