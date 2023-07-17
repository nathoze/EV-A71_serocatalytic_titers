
truncated_poisson_pmf <- function(lambda, truncation.point, x) {
  p <- dpois(x, lambda) /( ppois(truncation.point, lambda) -exp(-lambda))
  p[x > truncation.point] <- 0 # Set probability to 0 for values above truncation point
  p[x == 0] <- 0

  return(p)
}

decay_poisson_pmf <- function(lambda, truncation.point, x) {
  p <- dpois(x, lambda) /( ppois(truncation.point, lambda))
  p[x > truncation.point] <- 0 # Set probability to 0 for values above truncation point

  return(rev(p))
}

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

infection_probability <- function(foi,N.titers = 10){
  return(c(rep(1-exp(-foi),N.titers-1), 0))
}

titer_distribution <- function(transition.matrices){ # replace on specific indices

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
      titers = t(transition.matrices[[i+j]])%*% titers
      Titers.year[,i] = titers
    }
    j=j+1

    obs = get_observed_titers(Titers.year)

    df=rbind(df, data.frame( birth.year = birth.year,
                             age = rep(seq(1,n), each =titer.observable.max),
                             obs.proportion = as.numeric(obs), titer.class= rep(seq(1,titer.observable.max), n)))
  }

  df2 = df %>%
    mutate(sampling.year = age+birth.year) %>%
    #filter(age <= age.max) %>%
    filter(sampling.year >= min.year.sampling) %>%
    filter(sampling.year <= max.year.sampling)

  return(df2)

}

#
# titer_distribution_optimized <- function(all.parameters, year.changes = NULL){
#
#
#   # indices.changes = 1:N.FOI){ # replace on specific indices
#   birth.years = seq(min.year.sampling-age.max,max.year.sampling )
#   Y = birth.years
#   if(!is.null(year.changes)){
#     possible.foi = seq(birth.years[year.changes]-age.max, birth.years[year.changes])
#     Y = birth.years[!is.na(match(birth.years, possible.foi))]
#   }
#
#   df = NULL
#   j=0
#   for(birth.year in Y){
#
#     titers = Titers.0
#
#     years <- head(seq(birth.year,min(birth.year+age.max,max.year.sampling)),age.max)
#     n = length(years)
#     Titers.year = matrix(data = NA, nrow = N.titers, ncol = n)
#     i=0
#     for(foi.index in years){
#       i=i+1
#       titers = t(all.parameters$transition.matrices[[i+j]])%*% titers
#       Titers.year[,i] = titers
#     }
#     j=j+1
#
#     obs = get_observed_titers(Titers.year)
#
#     df=rbind(df, data.frame( birth.year = birth.year,
#                              age = rep(seq(1,n), each =titer.observable.max),
#                              obs.proportion = as.numeric(obs), titer.class= rep(seq(1,titer.observable.max), n)))
#   }
#
#   df2 = df %>%
#     mutate(sampling.year = age+birth.year) %>%
#     #filter(age <= age.max) %>%
#     filter(sampling.year >= min.year.sampling) %>%
#     filter(sampling.year <= max.year.sampling)
#
#   return(df2)
#
# }

get_all_parameters <- function(params){

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
update_all_parameters <- function(old.all.parameters, new.param, updated_index){

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
    new.all.parameters = get_all_parameters(new.all.parameters$params)
  }
  new.all.parameters$titer.distribution = titer_distribution(new.all.parameters$transition.matrices)

  return(new.all.parameters)

}

# Right-censor the data
get_observed_titers <- function(Titers.year, titer.observable.max = 5){

  M=as.matrix(Titers.year)
  if(ncol(M)==1){
    M[titer.observable.max] = sum(Titers.year[seq(titer.observable.max,nrow(Titers.year))])
    M = M[seq(1,titer.observable.max)]
  }else{
    M[titer.observable.max,] = colSums(Titers.year[seq(titer.observable.max,nrow(Titers.year)),])
    M = M[seq(1,titer.observable.max),]
  }

  return(M)
}

