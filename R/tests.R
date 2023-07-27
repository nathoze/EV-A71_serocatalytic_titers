params0  = c(  runif(n  = round(N.FOI/5), max = 0.7), runif(n = 1, min=0,max=5),1)
model_five_years$compute_loglik(all.params = model_five_years$get_all_parameters(params0))


P=model_constant$get_all_parameters(res$params[100,])


all.params <-  model_constant$update_all_parameters(P,new.param = res$params[100,1],updated_index =1 )

# model_constant$get_all_parameters(res$params[100,])

model_constant$compute_loglik(all.params = model_constant$get_all_parameters(params0))


set.seed(1)
inds_to_update=1
model_constant = list(compute_loglik = compute_loglik,
                      params0 = params0,
                      inds_to_update = inds_to_update,
                      is_invalid = is_invalid_model_constant,
                      get_all_parameters = get_all_parameters_model_constant,
                      update_all_parameters = update_all_parameters_model_constant)


res <-  run_MCMC_specify_model(model = model_constant,
                               mcmc_steps = 10,
                               mcmc_adaptive_steps = mcmc_adaptive_steps,
                               verbose = TRUE)
LL=c()
for(i in 1:10){
  LL[i] =  res$model$compute_loglik(res$model$get_all_parameters(res$params[i,]))
}
plot(LL)
lines(res$loglik)



