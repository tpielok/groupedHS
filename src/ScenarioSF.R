# replicate simulation study of Scheipl (2010)
# data generating process for medium sized AMs:
dgp_medium_am <- function(scenarios = list(),
                          dgp = 1,
                          isTestData = F) {
  p <- ifelse(scenarios$sparse[dgp] == "high", 20, 16)
  n <-  if (!isTestData) {
    scenarios$n[dgp]
  } else
    5e3
  snr <- scenarios$snr[dgp]
  corr <- scenarios$corr[dgp]
  
  
  # data
  x <- {
    S <- corr ^ as.matrix(dist(1:p))
    matrix(runif(n * p, -2, 2), nr = n) %*% chol(S)
  }
  colnames(x) <- paste0("x", 1:p)
  
  # true lin pred
  ff1 <- function (x)
    x
  ff2 <- function (x)
    x + (2 * x - 2) ^ 2 / 5.5
  ff3 <- function (x)
    - x + pi * sin(pi * x)
  ff4 <- function (x)
    .5 * x + 15 * (dnorm((x - .2) / .5) - dnorm(x + .4))
  
  f1 <- ff1(x[, 1])
  f2 <- ff2(x[, 2])
  f3 <- ff3(x[, 3])
  f4 <- ff4(x[, 4])
  
  if (scenarios$sparse[dgp] == "high") {
    f <- data.frame(cbind(f = f1 + f2 + f3 + f4, f1, f2, f3, f4))
    truthSm <- truthSmBoost <- c(F, T, T, T, rep(F, p - 4))
    truthFx <- truthFxBoost <- c(T, T, T, T, rep(F, p - 4))
    truth <-
      truthBoost <- c(T, F, T, T, T, T, T, T, rep(F, 2 * (p - 4)))
  } else {
    f5 <- 1.5 * ff1(x[, 5])
    f6 <- 1.5 * ff2(x[, 6])
    f7 <- 1.5 * ff3(x[, 7])
    f8 <- 1.5 * ff4(x[, 8])
    f9 <- 2 * ff1(x[, 9])
    f10 <- 2 * ff2(x[, 10])
    f11 <- 2 * ff3(x[, 11])
    f12 <- 2 * ff4(x[, 12])
    
    f <- cbind(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12)
    f <- data.frame(f)
    f$f <- rowSums(f)
    
    truthSm <- c(rep(c(F, T, T, T), t = 3), rep(F, p - 12))
    truthFx <- c(rep(T, 12), rep(F, p - 12))
    truth <-
      c(rep(c(T, F, T, T, T, T, T, T), t = 3), rep(F, 2 * (p - 12)))
  }
  
  
  formula <-
    formula(paste("y ~ ", paste(colnames(x), collapse = "+"), sep = ""))
  boost.formula <- formula(paste(
    "y ~ ",
    paste(
      "bols(",
      colnames(x),
      ")+bbs(",
      colnames(x),
      ",center=T, df=1, knots=10)",
      collapse = "+",
      sep = ""
    )
  ))
  true.gam.formula <- formula(paste("y ~ ",
                                    gsub(
                                      "(\\+){2,}",
                                      "",
                                      paste(
                                        ifelse(truthSm, "s(", ""),
                                        ifelse(truthSm |
                                                 truthFx, colnames(x), ""),
                                        ifelse(truthSm, ", k=12)", ""),
                                        collapse = "+",
                                        sep = ""
                                      )
                                    )))
  
  # response & test data
  eps <- scale(rnorm(n))
  eps <- sqrt(var(f$f) / snr) * eps
  y <- f$f + eps
  y_scal = 0
  
  test <- if (!isTestData) {
    dgp_medium_am(scenarios, dgp, isTestData = TRUE)
  } else
    NULL
  if (!isTestData) {
    # make sure learn and test data are rescaled the same way
    y <- scale(y)
    y_scal <- attr(y, "scaled:scale")
    f$f <-
      (f$f - attr(y, "scaled:center")) / attr(y, "scaled:scale")
    x <- scale(x) / 2
    test$f <-
      (test$f - attr(y, "scaled:center")) / attr(y, "scaled:scale")
    test$y <-
      (test$y - attr(y, "scaled:center")) / attr(y, "scaled:scale")
    test[, grep("x", colnames(test))] <-
      t((t(test[, grep("x", colnames(test))]) - attr(x, "scaled:center")) / (2 *
                                                                               attr(x, "scaled:scale")))
  }
  data <- data.frame(cbind(y, f, x))
  
  return(
    structure(
      data,
      test = test,
      formula = formula,
      boost.formula = boost.formula,
      true.gam.formula = true.gam.formula,
      y_scal = y_scal
    )
  )
}

getF <- function(data, num_cov, num_nz) {
  F = matrix(0, ncol = num_cov, nrow = nrow(data))
  for (i in 1:num_nz) {
    F[, i] = data[[paste0("f", as.character(i))]]
  }
  return(F)
}

createDatasets <- function(scenarios, reps){
  datasets = vector("list", nrow(scenarios))
  
  for (dgp in 1:nrow(scenarios)) {
    datasets[[dgp]] = vector("list", reps)
    for (j in 1:reps) {
      datasets[[dgp]][[j]] <-
        dgp_medium_am(scenarios = scenarios, dgp = dgp)
    }
  }
  return(datasets)
}


compute_ghs_mse <- function(num_cov, num_nz, data){
  psplines = vector("list", num_cov)
  knots = seq(-2, 2, length.out = 20)
  degree = 3
  d = 2
  
  
  for (i in 1:num_cov) {
    psplines[[i]] =  pspline1D(data[[paste0("x", as.character(i))]],
                               knots = knots,
                               degree = degree,
                               d = d)
  }
  
  
  sp = stan_params_sf(psplines, data$y, 1)
  if (scenarios$p[dgp] == "strict") {
    p0 = 1
  } else if (scenarios$p[dgp] == "opt") {
    p0 = num_nz * sp$num_params / sp$num_splines
  } else{
    p0 = sp$num_params - 1
  }
  
  tau0 = p0 / (sqrt(sp$num_data) * (sp$num_params - p0))
  
  sp = stan_params(psplines = psplines, Y=data$y, scale_global = tau0)

  ptm <- proc.time()
  fit <- rstan::sampling(attr(sp,"sm"),
                  data = sp,
                  iter = 500,
                  seed = cnt)
  elapsedTime = (proc.time() - ptm)[3]
  
  testdata = attr(data, "test")
  X = testdata[, grep("x", colnames(testdata))]
  
  esp = evalSP_sf(
    psplines = psplines,
    sp = sp,
    fit = fit,
    X = X
  )
  
  hs_mse = mean((esp$y_hat - testdata$y) ^ 2)
  attr(hs_mse, "gamma") = esp$gamma
  attr(hs_mse, "beta0") = esp$beta0
  attr(hs_mse, "f_hat") = esp$f 
  attr(hs_mse, "elapsedTime") = elapsedTime
  return(hs_mse)
}

compute_fo_mse <- function(f_hat, num_cov, num_nz, data){
  testdata = attr(data, "test")
  Fs = centerF(getF(testdata, num_cov, num_nz)) / attr(data, "y_scal")
  return(cmpEval(f_hat, Fs))
}

compute_ss_mse <- function (data){
  m <- spikeSlabGAM(attr(data, "formula"), data)
  ss_mse = mean((predict(m, attr(data, "test")) - attr(data, "test")$y) ^
                  2)
  return(ss_mse)
}

compute_gl_mse <- function(cnt, data){
  set.seed(cnt)
  X = data[, grep("x", colnames(data))]
  test = attr(data, "test")
  bases = pseudo.bases(X, degree = 10, df = 6)
  
  cvg =  cv.gamsel(X, data$y, bases = bases)
  
  X = test[, grep("x", colnames(test))]
  gl_mse = mean((test$y - predict(cvg$gamsel.fit, X)[, cvg$index.min]) ^
                                        2)
  return(gl_mse)
}

compute_oracle_mse <- function(data){
  # orakel-modell als benchmark:
  m_mgcv <- mgcv::gam(attr(data, "true.gam.formula"), data = data)
  oracle_mse = mean((predict(m_mgcv, newdata = attr(data, "test")) - attr(data, "test")$y) ^
                      2)
  return(oracle_mse)
}

if (TRUE
    ) {
  source("util-sf.R")
  
  library("spikeSlabGAM")
  library("rstan")
  library("gamsel")
  library("matlib")
  
  source("ghs-am.R")
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
    
  seed_cnt = 1
  set.seed(seed_cnt)
  dgp = 1
  
  resfile = get_sf_file()
  scenarios = get_sf_scenarios()
  
  reps = 100
  datasets = createDatasets(scenarios, reps)
  
  nsz =  nrow(scenarios)
  total = nrow(scenarios) * reps
  results = vector("list", 0)
  
  for (cnt in 1:total) {
    j = cnt %/% nsz + 1
    dgp   = cnt %% nsz + 1
    
    data = datasets[[dgp]][[j]]
    
    seed_cnt = dgp * reps + j
    set.seed(seed_cnt)
    
    num_cov <- ifelse(scenarios$sparse[dgp] == "high", 20, 16)
    num_nz <- ifelse(scenarios$sparse[dgp] == "high", 4, 12)
    
    hs_mse = compute_ghs_mse(num_cov, num_nz, data)
    quality = compute_fo_mse(attr(hs_mse,'f_hat'), num_cov, num_nz, data)
    
    result = list(
      quality = quality,
      j = j,
      dgp = dgp,
      hs_mse = c(hs_mse),
      gl_mse = compute_gl_mse(cnt,data),
      ss_mse = compute_ss_mse(data),
      oracle_mse = compute_oracle_mse(data),
      gamma = attr(hs_mse,'gamma'),
      beta0 = attr(hs_mse,'beta0'),
      elapsedTime = attr(hs_mse,'elapsedTime')
    )
    print(paste("Progress: ", as.character(cnt / total)))
    results[[as.character(cnt)]] = result
    save(results, file = resfile)
    
  }
}