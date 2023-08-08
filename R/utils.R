
run_MCMC <- function(model,
                     data,
                     proposal_type = NULL,
                     mcmc_steps = 1000, mcmc_adaptive_steps = 100,
                     sd_proposal = NULL, opt_acceptance = 0.24,
                     verbose = FALSE) {

  compute_loglik = model$compute_loglik
  is_invalid = model$is_invalid
  params0 = model$params0
  inds_to_update = model$inds_to_update
  update_all_parameters = model$update_all_parameters
  get_all_parameters = model$get_all_parameters
  fct_model_antibody_increase = model$fct_model_antibody_increase
  fct_model_antibody_decrease = model$fct_model_antibody_decrease
  n_params <- length(params0)

  # Preallocate containers
  loglik <- rep(0, mcmc_steps)
  accept <- matrix(0, ncol = n_params, nrow = mcmc_steps) # Acceptance rates
  params <- matrix(0, ncol = n_params, nrow = mcmc_steps) # Parameters
  params[1, ] <- params0

  print(get_all_parameters)
  all.params = get_all_parameters(params0, fct_model_antibody_increase,fct_model_antibody_decrease)

  old.all.params = all.params

  # Default choices when argument is NULL
  if (is.null(proposal_type)) { # Lognormal proposal for all parameters
    proposal_type = rep("lognorm", n_params)
  }
  if (is.null(inds_to_update)) { # Update all parameters
    inds_to_update <- 1:n_params
  }
  if (is.null(sd_proposal)) { # Start with 0.01 for all parameters
    sd_proposal <- rep(0.01, n_params)
  }

  # MCMC loop start
  step <- 1
  loglik[step] <- compute_loglik(all.params,data)

  accept[step, inds_to_update] <- 1
  accepted <- rep(0, n_params) # Number of accepted moves
  accepted[inds_to_update] <- 1
  for (step in 2:mcmc_steps) {

    if (verbose) print(step)
    params[step, ] <- params[step - 1, ]
    loglik[step] <- loglik[step - 1]

    for (k in inds_to_update) { # Parameter update start
      #    print(paste0("param ", k))
      old_param <- params[step, k]
      # Propose new param
      if (proposal_type[k] == "norm") {
        # Normal proposal
        new_param <- old_param + sd_proposal[k] * rnorm(1)
        log_proposal <- 0.0
      } else {
        # Lognormal proposal (default)
        new_param <- old_param * exp(sd_proposal[k] * rnorm(1))
        log_proposal <- log(new_param) - log(old_param)
      }
      # Reject immediately if new param is invalid
      if (is_invalid(k, new_param)) {
        next
      }

      # Update the transformed parameters
      new.all.params <-  update_all_parameters(old.all.params,
                                               new.param = new_param,
                                               updated_index = k,
                                               fct_model_antibody_increase,
                                               fct_model_antibody_decrease)

      params[step, ] = new.all.params$params

      new_loglik <- compute_loglik(new.all.params,data)

      # Metropolis-Hastings
      log_ratio <- new_loglik - loglik[step] + log_proposal


      #print(log_ratio)
      if(!is.na(log_ratio)){
        if (log(runif(1)) < log_ratio) {
          loglik[step] <- new_loglik
          accepted[k] <- accepted[k] + 1
          old.all.params= new.all.params
        } else {
          new.all.params = old.all.params
          params[step,k] <- old_param
        }
      } else {
        new.all.params = old.all.params
        params[step,k] <- old_param

      }
    } # Parameter update end

    # Acceptance rates

    accept[step, ] <- accepted / step

    # Adjust proposal sd to get closer to 'opt_acceptance'
    if (step <= mcmc_adaptive_steps) {
      delta <- 1 - step / mcmc_adaptive_steps
      for (k in inds_to_update) {
        diff <- accept[step, k] - opt_acceptance
        sd_proposal[k] <- sd_proposal[k] * (1.0 + delta * diff)
        #sd_proposal[k] = 0.01
      }
    }

  } # MCMC loop end

  list(loglik = loglik, params = params, accept = accept, model = model, sd_proposal= sd_proposal, data = data)
}


define_model <- function(fct_model_antibody_increase = get_increase_matrix,
                         fct_model_antibody_decrease = get_decay_matrix,
                         model_foi =  'constant',
                         compute_loglik,
                         params0,
                         inds_to_update){

  source('R/FOI_models.R') # when sourcing we redefine the different functions

  if(model_foi == 'constant'){
    is_invalid = is_invalid_model_constant
    get_all_parameters = get_all_parameters_model_constant
    update_all_parameters = update_all_parameters_model_constant
  }
  if(model_foi == '5years'){
    is_invalid = is_invalid_model_five_years
    get_all_parameters = get_all_parameters_model_five_years
    update_all_parameters = update_all_parameters_model_five_years
  }
  if(model_foi == 'independent'){
    is_invalid = is_invalid_model_independent
    get_all_parameters = get_all_parameters_model_independent
    update_all_parameters = update_all_parameters_model_independent
  }
  if(model_foi == 'peak_constant'){
    is_invalid = is_invalid_model_peak_constant
    get_all_parameters = get_all_parameters_model_peak_constant
    update_all_parameters = update_all_parameters_model_peak_constant
  }

  model = list(fct_model_antibody_increase = fct_model_antibody_increase,
               fct_model_antibody_decrease = fct_model_antibody_decrease,
               compute_loglik = compute_loglik,
               params0 = params0,
               inds_to_update = inds_to_update,
               is_invalid = is_invalid,
               get_all_parameters = get_all_parameters,
               update_all_parameters = update_all_parameters)
  return(model)
}

compute_DIC <- function(results, burn_in){

  #https://en.wikipedia.org/wiki/Deviance_information_criterion

  data = results$data
  params = results$params[-seq(1,burn_in), ]

  avg_params= colMeans(params)
  all.params = results$model$get_all_parameters(params = avg_params,
                                                fct_model_antibody_increase = results$model$fct_model_antibody_increase,
                                                fct_model_antibody_decrease = results$model$fct_model_antibody_decrease)

  loglik = results$loglik[-seq(1,burn_in)]

  avg_deviance = mean(-2*loglik)
  deviance_avg_params = -2*results$model$compute_loglik(all.params, data)
  pD = avg_deviance-deviance_avg_params
  DIC = pD + avg_deviance
  return(list(pD = pD,
              DIC = DIC))
}


quantile025 <- function(X){
  return(as.numeric(quantile(X, probs=0.025 )))
}

quantile975 <- function(X){
  return(as.numeric(quantile(X, probs=0.975 )))
}


# Log likelihood
compute_loglik <- function(all.params,data) {

  prob.titers =  left_join(all.params$titer.distribution,
                           data,
                           by = c("birth.year", "age", "titer.class", "sampling.year"))  %>%
    group_by(age, sampling.year) %>%
    summarise(ll  = dmultinom(x = n,prob = obs.proportion, log=TRUE), .groups='drop')  %>%
    summarise(s= sum(ll))

  return(prob.titers$s)

}
