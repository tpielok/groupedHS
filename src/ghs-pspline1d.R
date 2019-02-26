library("splines")
library("matlib")

diffMatrix <- function(d, n) {
  D = matrix(0, n - d, n)
  for (i in 0:d) {
    D = D + diag(n + i)[(i + 1):(n - d + i), 1:n] *
      choose(d, i) * (-1) ^ ((d + 1) %/% 2 + i)
  }
  return(D)
}

penaltyMatrix <- function(d, n) {
  #todo vector edge case
  D = diffMatrix(d, n)
  return (t(D) %*% D)
}

range01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}


transformInput <- function(x,trafo){
  return(t((t(x) - trafo[2,])/trafo[1,]))
}

pspline1D <- function(x, knots = NULL, degree = 3, d = 2) {
  if(is.null(knots)){
    knots = seq(0, 1, length.out = 20)
    trafo = c(max(x)-min(x),min(x))
    x = transformInput(x,matrix(trafo))
  }else{
    trafo = c(1,0)
  }
  
  B = bs(x,
         knots = knots,
         degree = degree,
         intercept = FALSE)
  
  p = ncol(B)
  iP = diag(p)
  iP[(1+d):p,] = diffMatrix(d,p) 
  P = inv(iP) # parametrization matrix
  B = B %*% P # transform desing matrix
  
  return(list(
    B = B,
    P = P,
    knots = knots,
    degree = degree,
    trafo = trafo
  ))
}

collapseMatrices <- function(psplines, z) {
  num = length(psplines)
  Offset = 1
  Start = vector("numeric", num)
  End   = vector("numeric", num)
  collapsed = vector("numeric", 0)
  
  for (i in 1:num) {
    pspline = psplines[[i]]
    Start[i] = Offset
    Offset = Offset + length(psplines[[1]][[z]])
    End[i] = Offset - 1
    collapsed = c(collapsed, as.vector(pspline[[z]]))
  }
  return(list(
    Start = matrix(Start, nrow = 1),
    End = matrix(End, nrow = 1),
    collapsed = collapsed
  ))
}

createParams <- function(psplines) {
  num = length(psplines)
  Offset = 1
  Start = vector("numeric", num)
  End   = vector("numeric", num)
  Len   = vector("numeric", num)
  
  for (i in 1:num) {
    pspline = psplines[[i]]
    Start[i] = Offset
    Len[i] = ncol(psplines[[i]][['B']])
    Offset = Offset + Len[i]
    End[i] = Offset - 1
  }
  return(list(
    Start = matrix(Start, nrow = 1),
    End = matrix(End, nrow = 1),
    Len = matrix(Len, nrow = 1)
  ))
}

stan_params_sf <- function(psplines, Y, scale_global, nu_local=1, nu_global=1) {
  cBs = collapseMatrices(psplines, 'B')
  gamma = createParams(psplines)
  trafo = sapply(psplines, function(x) x$trafo)
  
  return(
    list(
      num_data = length(Y),
      Y = Y,
      num_splines = length(psplines),
      num_params = sum(gamma$Len),
      len_Bs = length(cBs$collapsed),
      bStart = cBs$Start,
      bEnd = cBs$End,
      Bs = cBs$collapsed,
      gStart = gamma$Start,
      gLen = gamma$Len,
      gEnd = gamma$End,
      nu_global = nu_global,
      nu_local = nu_local,
      scale_global = scale_global,
      trafo = trafo
    )
  )
}

evalSP_sf <- function(psplines, sp, fit, X) {
  la = rstan::extract(fit, permuted = TRUE)
  
  u = colMeans(la$u)
  beta0 = mean(la$beta0)
  gammas = vector("numeric",length = 0)
  f = matrix(nrow = nrow(X), ncol = sp$num_splines)
  
  transformInput(X, sp$trafo)
  
  for (i in 1:length(psplines)) {
    p = psplines[[i]]
    B = bs(
      X[, i],
      knots = p$knots,
      degree = p$degree,
      intercept = FALSE
    )
    gamma = p$P %*% u[sp$gStart[1, i]:sp$gEnd[1, i]]
    f[, i] = B %*% gamma
    mf = mean(f[, i])
    f[, i] = f[, i] - mf
    beta0 = beta0 + mf
    gammas = c(gammas,gamma)
  }
  return(list(
    gamma = gammas,
    f = f,
    beta0 = beta0,
    y_hat = (rowSums(f) + beta0)
  ))
}

centerF <- function(f_true){
  return(t(t(f_true) - colMeans(f_true)))
}

cmpEval <- function (f_hat, f_true){
  return(colMeans((f_hat - centerF(f_true))^2))
}

if (FALSE) { # Example
  library("rstan")
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  sm = stan_model("ghs-pspline1d.stan")
  
  cnt_seed = 1
  
  set.seed(cnt_seed)
  num_data = 200
  num_cov = 3
  X <-
    matrix(runif(
      n = num_data * num_cov,
      min = -5,
      max = 5
    ), ncol = num_cov)
  
  f_true = matrix(0,nrow=num_data, ncol=1)  
  f_true = matrix(0,nrow=num_data, ncol=3)
  f_true[,1] = sin(X[, 1])
  f_true[,2] = cos(X[, 2])
  
  Y_true <- f_true[,1] + f_true[,2]
  
  Y <- Y_true + rnorm(length(Y_true), 0, .2) # adding noise
  
  degree = 3
  d = 3
  knots = seq(-5, 5, length.out = 20)
  p1 <- pspline1D(X[, 1],
                  knots = knots,
                  degree = degree,
                  d = d)
  p2 <- pspline1D(X[, 2],
                  knots = knots,
                  degree = degree,
                  d = d)
  p3 <- pspline1D(X[, 3],
                  knots = knots,
                  degree = degree,
                  d = d)
  
  
  psplines = list(p1, p2, p3)
  sp = stan_params_sf(psplines, Y)
  
  fit = sampling(
    sm,
    data = sp,
    iter = 500,
    control = list(adapt_delta = 0.999),
    seed = cnt_seed,
    sample_file = tempdir()
  )
  
  la = extract(fit, permuted = TRUE)
  plot(colMeans(la$Y_hat) - Y_true)
  
  X <-
    matrix(runif(
      n = num_data * num_cov,
      min = -5,
      max = 5
    ), ncol = num_cov)
  
  f_true = matrix(0,nrow=num_data, ncol=1)  
  f_true = matrix(0,nrow=num_data, ncol=3)
  f_true[,1] = sin(X[, 1])
  f_true[,2] = cos(X[, 2])
  
  Y_true <- f_true[,1] + f_true[,2]  
  
  esp = evalSP_sf(
    psplines = psplines,
    sp = sp,
    fit = fit,
    X = X
  )
  
  plot(esp$y_hat - Y_true)
  mean((esp$y_hat - Y_true) ^ 2)
  mean((esp$y_hat - Y)^2)
  
  cf_true = centerF(f_true)
  plot(X[, 1], esp$f[, 1])
  points(X[, 1], cf_true[,1] , col = 'blue') 
  plot(X[, 2], esp$f[, 2])
  points(X[, 2], cf_true[,2] , col = 'blue') 
  plot(X[, 3], esp$f[, 3])
  points(X[, 3], cf_true[,3] , col = 'blue') 
  
  cmpEval(esp, f_true)
  
  plot(density(rowSums(esp$f) + esp$beta0 - Y_true))
}