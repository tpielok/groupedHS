glevels = c(3, 5, 9)         # how many levels per group (small, medium, large)
num_dg  = length(glevels)    # number of different groups
beta_min = 1
beta_max = 2

createGroupedCoefs <- function(dgp, cur_sp_level) {
  # create betas depending on group size and level of sparsity
  templ_betas = vector("numeric", 0)
  templ_group_ids = vector("numeric", 0)
  
  i = 1
  for (gl in glevels) {
    templ_betas = c(templ_betas, seq(beta_min, beta_max, length.out = gl - 1))
    templ_group_ids = c(templ_group_ids, rep(i, gl - 1))
    i = i + 1
  }
  
  num_nonzero = length(templ_betas)
  betas = vector("numeric", 0)
  group_ids = vector("numeric", 0)
  
  betas = c(templ_betas, rep(0, num_nonzero * cur_sp_level))
  
  group_ids = templ_group_ids
  for (i in 1:cur_sp_level) {
    group_ids = c(group_ids, templ_group_ids + i * num_dg)
  }
  
  return(list(
    betas = betas,
    group_ids = group_ids,
    num_nonzero = num_nonzero
  ))
}

dgp_lp <- function(scenarios,
                   dgp = 1,
                   seed_cnt,
                   isTestData = F) {
  p <- ifelse(scenarios$sparse[dgp] == "high", 20, 16)
  n <-  if (!isTestData) {
    scenarios$n[dgp]
  } else
    5e3
  
  cur_sp_level = ifelse(scenarios$sparse[dgp] == "high", 2, 1)
  gCoef = createGroupedCoefs(dgp, cur_sp_level)
  D = length(gCoef$betas)
  
  p0 = c(1, gCoef$num_nonzero, D - 1)
  num_groups = (cur_sp_level + 1) * num_dg
  
  if (scenarios$p[dgp] == "strict") {
    cur_p0 = 1
  } else if (scenarios$p[dgp] == "opt") {
    cur_p0 = gCoef$num_nonzero
  } else if (scenarios$p[dgp] == "loose") {
    cur_p0 = D - 1
  }
  
  # data
  set.seed(seed_cnt)
  inpMat = matrix(nrow = n, ncol = num_dg * (1 + cur_sp_level))
  for (i in 1:ncol(inpMat)) {
    cur_glevel = glevels[((i - 1) %% num_dg) + 1]
    # todo: change to draw till valid
    inpMat[, i] = as.character(c(
      sample(1:cur_glevel, cur_glevel, FALSE),
      sample(1:cur_glevel, n - cur_glevel, TRUE)
    ))
  }
  
  df = data.frame(inpMat)
  mM = model.matrix(~ ., df)
  
  x = mM[, 2:ncol(mM)] # intercept gets estimated seperately
  y = x %*% gCoef$betas
  
  tau0 = cur_p0 / (sqrt(n) * (D - cur_p0))
  
  # response & test data
  # adjust snr
  eps = scale(rnorm(n))
  eps = sqrt(var(y)[1, 1] / scenarios$snr[dgp]) * eps
  y = data.frame(y = y + eps)
  
  test <- if (!isTestData) {
    dgp_lp(scenarios,
           dgp,
           seed_cnt = seed_cnt + 1,
           isTestData = TRUE)
  } else
    NULL
  
  data <- data.frame(cbind(y, x))
  
  return(
    structure(
      data,
      test = test,
      tau0 = tau0,
      group_ids = gCoef$group_ids,
      num_groups = num_groups,
      df_tp = df[, 1:num_dg],
      # true predictors
      beta = gCoef$betas,
      num_nonzero = gCoef$num_nonzero
    )
  )
}

createDatasets <- function(reps, scenarios) {
  datasets = vector("list", nrow(scenarios))
  i = 1
  for (dgp in 1:nrow(scenarios)) {
    datasets[[dgp]] = vector("list", reps)
    for (j in 1:reps) {
      datasets[[dgp]][[j]] <-
        dgp_lp(scenarios = scenarios, dgp = dgp, i)
      i = i + 1
    }
  }
  return(datasets)
}

compute_mse <- function(ds) {
  sp = stan_params(
    Y = ds$y,
    X = ds[, -which(names(ds) %in% c("y"))],
    group_ids = attr(ds, "group_ids"),
    num_groups = attr(ds, "num_groups"),
    scale_global = attr(ds, "tau0")
  )
  
  ptm <- proc.time()
  fit = sampling(
    attr(sp, 'sm'),
    data = sp,
    iter = 500,
    control = list(adapt_delta = 0.99),
    seed = seed_cnt,
    sample_file = tempdir()
  )
  elapsedTime = (proc.time() - ptm)[3]
  
  test_data = attr(ds, "test")
  ev = evalSP_lp(sp = sp,
                 fit = fit,
                 X = as.matrix(test_data[, -which(names(ds) %in% c("y"))]))
  ghs_mse = mean((ev$y_hat - test_data$y) ^ 2)
  attr(ghs_mse, 'elapsedTime') = elapsedTime
  attr(ghs_mse, 'beta') = ev$beta
  attr(ghs_mse, 'beta0') = ev$beta0
  return(ghs_mse)
}

compute_param_mse <- function(ds, beta) {
  num_nonzero = attr(ds, "num_nonzero")
  cur_sp_level = ifelse(scenarios$sparse[dgp] == "high", 2, 1)
  sq_res = (beta - attr(ds, "beta")) ^ 2
  
  mse_nonzero = mean(sq_res[1:num_nonzero])
  mse_zero    = mean(sq_res[(num_nonzero + 1):length(sq_res)])
  
  mse_levels = matrix(nrow = 2, ncol = num_dg) # 1: nonzero 2: zero
  
  sq_res = matrix(sq_res, nrow = cur_sp_level + 1, byrow = TRUE)
  offset = 1
  for (l in 1:num_dg) {
    mse_levels[1, l] = sqrt(mean(sq_res[1, offset:(offset + glevels[l] - 2)]))
    
    mse_levels[2, l] =  sqrt(mean(colMeans(matrix(sq_res[2:(cur_sp_level + 1), offset:(offset + glevels[l] - 2)], nrow =
                                                    cur_sp_level))))
    offset = offset + glevels[l] - 1
  }
  return(list(
    mse_nonzero = mse_nonzero,
    mse_zero = mse_zero,
    mse_levels = mse_levels
  ))
}

compute_oracle_mse <- function(ds) {
  test_data = attr(ds, "test")
  m = lm(ds$y ~ ., attr(ds, "df_tp"))
  oracle_y = predict(m, attr(test_data, "df_tp"))
  oracle_mse = mean((test_data$y - oracle_y) ^ 2)
}

if (TRUE) {
  source("ghs-am.R")
  source("util-lp.R")
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  reps = 1
  seed_cnt = 1
  set.seed(seed_cnt)
  
  scenarios = get_lp_scenarios()
  resfile = get_lp_file()
  datasets = createDatasets(reps, scenarios)
  
  nsz =  nrow(scenarios)
  total = nsz * reps
  
  results = vector("list", 0)
  
  for (cnt in 1:total) {
    j = cnt %/% nsz + 1       # number of repetition
    dgp   = cnt %% nsz + 1
    
    ds = datasets[[dgp]][[j]]
    
    ghs_mse = compute_mse(ds)
    pmses = compute_param_mse(ds, attr(ghs_mse, 'beta'))
    oracle_mse = compute_oracle_mse(ds)
    
    result = list(
      mse_levels = pmses$mse_levels,
      mse_zero = pmses$mse_zero,
      mse_nonzero = pmses$mse_nonzero,
      j = j,
      dgp = dgp,
      ghs_mse = ghs_mse,
      oracle_mse = oracle_mse,
      beta = attr(ghs_mse, 'beta'),
      beta0 = attr(ghs_mse, 'beta0'),
      elapsedTime = attr(ghs_mse, 'elapsedTime')
    )
    print(paste("Progress: ", as.character(cnt / total)))
    results[[as.character(cnt)]] = result
    save(results, file = resfile)
  }
}
