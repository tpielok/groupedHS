library("splines")
library("rstan")
library("taRifx")
source("ghs-pspline1d.R")
source("ghs-lp.R")

if(!exists("sm_lp") || !exists("sm_sf") || !exists("sm_am")){
  sm_lp =  stan_model("ghs-lp.stan")
  sm_sf =  stan_model("ghs-pspline1d.stan")
  sm_am =  sm = stan_model("ghs-am.stan")
}

stan_params <-
  function(Y,
           psplines = NULL,
           X = NULL,
           group_ids = c(0),
           num_groups = 0,
           scale_global = 1,
           nu_local = 1,
           nu_global = 1) {
    num_data <- length(Y)
    
    if (!is.null(psplines)) {
      sp_sf = stan_params_sf(
        psplines,
        Y,
        scale_global = scale_global,
        nu_local = nu_local,
        nu_global = nu_global
      )
    }
    if (!is.null(X)) {
      sp_lp = stan_params_lp(
        Y,
        X,
        group_ids,
        num_groups,
        scale_global = scale_global,
        nu_local = nu_local,
        nu_global = nu_global
      )
    }
    
    if(!is.null(psplines) && !is.null(X)){
      sp = merge.list(sp_sf,sp_lp)
      sm = sm_am
    }else if(!is.null(psplines)){
      sp = sp_sf
      sm = sm_sf
    } else if(!is.null(X)){
      sp = sp_lp
      sm = sm_lp
    }
    attr(sp, 'sm') = sm
    return(sp)
}

if (FALSE){ # Example
  n = 100
  x1 = runif(n)
  x2 = sample(LETTERS[1:3], replace = T, size = n)
  y = sin(4*x1) + ifelse(x2 == "B", 0.5,0) + ifelse(x2 == "C", -0.5,0) 
  
  df = data.frame(x1,x2,y)
  
  mm = model.matrix(y ~ x2, df)
  X = mm[,2:3]
  group_ids = attr(mm, "assign")[2:3]
  num_groups = length(unique(group_ids))
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  degree = 3
  d = 3
  knots = seq(0, 1, length.out = 20)
  p1 <- pspline1D(df$x1,
                  knots = knots,
                  degree = degree,
                  d = d)
  sp = stan_params(Y=df$y,psplines = list(p1),X = X,group_ids = group_ids,num_groups = num_groups,
                   scale_global = 1)
  
  fit = sampling(
    attr(sp,'sm'),
    data = sp,
    iter = 200,
    control = list(adapt_delta = 0.999),
    seed = 1,
    sample_file = tempdir()
  )
  
  la = extract(fit, permuted = TRUE)
  plot(colMeans(la$Y_hat) - df$y)
  colMeans(la$beta)
  ev = evalSP_sf(list(p1),sp,fit,data.frame(df$x1))
  plot(df$x1,ev$y_hat)
}
