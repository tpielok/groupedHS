stan_params_lp <-
  function(Y,
           X,
           group_ids,
           num_groups,
           scale_global,
           nu_global = 1,
           nu_local = 1) {
    return(
      list(
        num_data = length(Y),
        Y = Y,
        X = X,
        num_groups = num_groups,
        num_linparam = ncol(X),
        group_ids = matrix(group_ids),
        scale_global = scale_global,
        nu_global = nu_global,
        nu_local = nu_local
      )
    )
    
  }

evalSP_lp <- function(sp, fit, X) {
  la = rstan::extract(fit, permuted = TRUE)
  
  beta0 = mean(la$beta0)
  beta = colMeans(la$beta)
  
  return(list(
    beta = beta,
    beta0 = beta0,
    y_hat = X %*% beta + beta0
  ))
}