source("ghs-am.R")

library("mlbench")
library("gamsel")
library("data.table")
library("MASS")
library("dismo")
library("plyr")
library("dplyr")
library("spikeSlabGAM")
library("rstan")

getBoston <- function(){
  tmp = new.env()
  data("BostonHousing", envir = tmp)
  return(tmp[["BostonHousing"]])
}

prepareBoston <- function(bn){
  ret = data.frame(bn %>% select(-c("chas","rad","medv")) %>% sapply(range01) )
  ret$medv = bn$medv
  return(ret)
}

createRepeatedCV <- function(bn, folds, repeats) {
  cv_design = replicate(repeats, kfold(bn, folds))
  return(list(
    bn = bn,
    cv_design = cv_design,
    repeats = repeats,
    folds = folds
  ))
}

ghs_mses <- function(xtrain, xtest, train, test){
  psplines = vector("list", length = ncol(train) - 1)
  i = 1
  degree = 3
  d = 2
  for (x in xtrain) {
    psplines[[i]] =  pspline1D(x,
                               knots = seq(0, 1, length.out = 20),
                               degree = degree,
                               d = d)
    i = i + 1
  }
  
  sp_g = stan_params(psplines = psplines, Y = train$medv, scale_global = 1)
  
  p0s = c(0.1 * sp_g$num_params,
          0.2 * sp_g$num_params,
          0.4 * sp_g$num_params)
  mses = vector("numeric", length = length(p0s))
  
  i = 1
  for (p0 in p0s) {
    tau0 = p0 / (sqrt(sp_g$num_data) * (sp_g$num_params - p0))
    sp = stan_params( psplines = psplines, Y= train$medv, scale_global = tau0)
    
    ptm <- proc.time()
    fit <- sampling(
      attr(sp,'sm'),
      data = sp,
      iter = 100,
      seed = floor(p0)
    )
    elapsedTime = (proc.time() - ptm)[3]
    
    esp = evalSP_sf(
      psplines = psplines,
      sp = sp,
      fit = fit,
      X = xtest
    )

    mses[i]  = mean((esp$y_hat - test$medv) ^ 2)
    i = i + 1
  }
  
  return(mses)
}

gamsel_mse <- function(xtrain, xtest, train, test){
  cvg =  cv.gamsel(xtrain, train$medv)
  gl_mse = mean((test$medv - predict(cvg$gamsel.fit, xtest)[, cvg$index.min]) ^
                  2)
  return(gl_mse)
}

ss_mse <- function(xtrain, xtest, train, test){
  formula <-
    formula(paste("medv ~ ", paste(colnames(xtrain), collapse = "+"), sep = ""))
  m <- spikeSlabGAM(formula, train)
  
  ss_mse = mean((predict(m, xtest) - test$medv) ^
                  2)
  return(ss_mse)
}

lm_mse <- function(train, test){
  m = lm(train, formula = medv ~ .)
  lm_mse = mean((test$medv - predict(m, test)) ^ 2)
  return(lm_mse)
}

gam_mse <- function(train, test, xtrain){
  gam.formula <- formula(paste("medv ~ ",
                               gsub(
                                 "(\\+){2,}",
                                 "",
                                 paste(
                                   "s(",
                                   colnames(xtrain),
                                   ", k=6)",
                                   collapse = "+",
                                   sep = ""
                                 )
                               )))
  m_mgcv <- mgcv::gam(gam.formula, data = train)
  gam_mse = mean((predict(m_mgcv, newdata = test) - test$medv) ^
                   2)
  return(gam_mse)
}

computeMSEs <- function(cv, k, bn){
  train = bn %>% filter(!(cv == k))
  test = bn %>% filter(cv == k)
  xtrain = train %>% select(-"medv")
  xtest = test %>% select(-"medv")
  
  return(list(
    ghs_mses = ghs_mses(xtrain = xtrain, xtest = xtest, train = train, test = test),
    gs_mse = gamsel_mse(xtrain = xtrain, xtest = xtest, train = train, test = test),
    ss_mse = ss_mse(xtrain = xtrain, xtest = xtest, train = train, test = test),
    lm_mse = lm_mse(train, test),
    gam_mse = gam_mse(train,test,xtrain)
  ))
}

doCV <- function(cv, folds, bn){
  results = sapply(seq(1,folds), function(x) computeMSEs(cv, x, bn))
  return(results)
}

doRepeatedCV <- function(rcv) {
  results = unlist(lapply(data.frame(rcv$cv_design), function(x)
    doCV(x, rcv$folds, rcv$bn)))
  
  cnames = c("p1","p2","p3","gs","ss","lm","gam")
  dfres = data.frame(matrix(as.vector(results),ncol=length(cnames),byrow = T))
  colnames(dfres) = cnames
  
  return(list(results = results,
              sresults = dfres))
}

permute <- function(x){
  sample(x,length(x),replace = F)
}

if(TRUE){
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  bn = getBoston() %>% prepareBoston()
  
  sparse_bn = bn %>% select(-"medv") %>% 
    mutate_if(is.numeric,funs('p' = permute(.)))
  sparse_bn$medv = bn$medv
  
  set.seed(1)
  rcv = createRepeatedCV(bn, folds = 5, repeats = 10)
  results = doRepeatedCV(rcv)
  save(results, file="Boston.RData")
  
  set.seed(2)
  rcv = createRepeatedCV(sparse_bn, folds = 5, repeats = 10)
  sparse_results = doRepeatedCV(rcv)
  save(sparse_results, file="sparseBoston.RData")
}
