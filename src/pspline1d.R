library("splines")

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

pspline1D <- function(x, knots, degree = 3, d = 2) {
  B = bs(x,
         knots = knots,
         degree = degree,
         intercept = FALSE)
  
  p = ncol(B)
  P = penaltyMatrix(d, p)
  rankP = sum(p - d)
  
  return(list(
    B = B,
    P = P,
    rankP = rankP,
    knots = knots,
    degree = degree
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

stan_params <- function(psplines, num_data, Y) {
  cBs = collapseMatrices(psplines, 'B')
  cPs = collapseMatrices(psplines, 'P')
  
  gamma = createParams(psplines)
  
  pRanks = vector("numeric", 0)
  for (pspline in psplines) {
    pRanks = c(pRanks, pspline$rankP)
  }
  
  return(
    list(
      num_data = num_data,
      Y = Y,
      num_splines = length(psplines),
      num_params = sum(gamma$Len),
      len_Bs = length(cBs$collapsed),
      len_Ps = length(cPs$collapsed),
      bStart = cBs$Start,
      bEnd = cBs$End,
      Bs = cBs$collapsed,
      pStart = cPs$Start,
      pEnd = cPs$End,
      Ps = cPs$collapsed,
      gStart = gamma$Start,
      gLen = gamma$Len,
      gEnd = gamma$End,
      pRanks = matrix(pRanks, nrow = 1),
      nu_global = 1,
      nu_local = 1,
      scale_global = 1
    )
  )
}

evalSP <- function(psplines, sp, fit, X) {
  la = extract(fit, permuted = TRUE)
  
  gamma = colMeans(la$gamma)
  beta0 = mean(la$beta0)
  
  f = matrix(nrow = sp$num_data, ncol = sp$num_splines)
  for (i in 1:length(psplines)) {
    p = psplines[[i]]
    B = bs(
      X[, i],
      knots = p$knots,
      degree = p$degree,
      intercept = FALSE
    )
    f[, i] = B %*% gamma[sp$gStart[1, i]:sp$gEnd[1, i]]
    mf = mean(f[, i])
    f[, i] = f[, i] - mf
    beta0 = beta0 + mf
  }
  return(list(
    f = f,
    beta0 = beta0,
    y_hat = (rowSums(f) + beta0)
  ))
}

centerF <- function(f_true){
  return(t(t(f_true) - colMeans(f_true)))
}

cmpEval <- function (esp, f_true){
  return(colMeans((esp$f - centerF(f_true))^2))
}

if (FALSE) {
  library("rstan")
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  sm = stan_model("pspline1d.stan")
  
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
  
  f_true = matrix(0,nrow=num_data, ncol=3)
  f_true[,1] = sin(X[, 1])
  f_true[,2] = cos(X[, 2])
  
  Y_true <- f_true[,1] + f_true[,2]
  
  Y <- Y_true + rnorm(length(Y_true), 0, .2) # adding noise
  
  degree = 2
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
  sp = stan_params(psplines, num_data, Y)
  
  fit = sampling(
    sm,
    data = sp,
    iter = 500,
    control = list(adapt_delta = 0.999),
    seed = cnt_seed,
    sample_file = tempdir()
  )
  
  X <-
    matrix(runif(
      n = num_data * num_cov,
      min = -5,
      max = 5
    ), ncol = num_cov)
  
  f_true = matrix(0,nrow=num_data, ncol=3)
  f_true[,1] = sin(X[, 1])
  f_true[,2] = cos(X[, 2])
  
  Y_true <- f_true[,1] + f_true[,2]  
  
  esp = evalSP(
    psplines = psplines,
    sp = sp,
    fit = fit,
    X = X
  )
  
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
