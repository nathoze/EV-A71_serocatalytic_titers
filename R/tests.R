params0  = c(  runif(n  = round(N.FOI/5), max = 0.7), runif(n = 1, min=0,max=5),1)
model_five_years$compute_loglik(all.params = model_five_years$get_all_parameters(params0))


P=model_constant$get_all_parameters(res$params[100,])


all.params <-  model_constant$update_all_parameters(P,new.param = res$params[100,1],updated_index =1 )

# model_constant$get_all_parameters(res$params[100,])

model_constant$compute_loglik(all.params = model_constant$get_all_parameters(params0))
