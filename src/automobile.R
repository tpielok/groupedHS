source("ghs-am.R")

library("data.table")
library("dismo")
library("plyr")
library("dplyr")
library("rstan")
library("gglasso")
library("spikeSlabGAM")

downloadAM <- function(auto_url) {
  automobile <-
    fread(auto_url,
          na.strings = "?")
  return(automobile)
}

prepareAM <- function(am) {
  # Only select categorical variables and numerical output
  sam = am %>% select(c(3:9, 15:16, 18, 26))
  colnames(sam) = c(
    "make",
    "fuel-type",
    "aspiration",
    "num-of-doors",
    "body-style",
    "drive-wheels",
    "engine-location",
    "engine-type",
    "num-of-cylinders",
    "fuel-system",
    "price"
  )
  
  pAM =  data.frame(select(sam,-"price") %>% sapply(function(x)
    as.factor(x)))
  pAM$price = sam$price
  return(pAM %>% filter(complete.cases(.)))
}

getCounts <- function(cnames, am) {
  return(lapply(X = cnames, function(x)
    count(am, am[, x]) %>% rename(!!x := "am[, x]")))
}

goodShuffleSplit  <- function(all_counts, test_counts) {
  for (i in 1:length(all_counts)) {
    ac = all_counts[[i]]
    mc = merge(ac, test_counts[[i]], by = colnames(ac)[1])
    if (min(mc$n.x - mc$n.y) <= 0) {
      return(F)
    }
  }
  return(T)
}

shuffleSplit <- function(df, ptrain) {
  return(c(rep(1, ceiling(ptrain * nrow(
    df
  ))), rep(2, floor((1 - ptrain) * nrow(df)
  ))) %>%
    sample(size = nrow(df), replace = F))
}

shuffleSplitLevel <- function(am, ptrain) {
  # Guarantee that that the test data does not contain unknown levels
  cnames =  colnames(select(am,-"price"))
  all_counts = getCounts(cnames, am)
  
  badSS = T
  while (badSS) {
    badSS = F
    ss = shuffleSplit(am, ptrain)
    sam = am %>% filter(ss == 2)
    if (!goodShuffleSplit(all_counts, getCounts(cnames, sam))) {
      badSS = T
    }
  }
  return(ss)
}


createShuffleSplits <- function(am, ptrain, repeats) {
  ss_design = replicate(repeats, shuffleSplitLevel(am, ptrain))
  return(list(
    am = am,
    ss_design = ss_design,
    repeats = repeats
  ))
}

ghs_mse <-
  function(X_train,
           X_test,
           y_train,
           y_test,
           group_ids,
           num_groups,
           tau0) {
    
    sp = stan_params(
      Y = y_train,
      X = X_train,
      group_ids = group_ids,
      num_groups = num_groups,
      scale_global = tau0
    )
    
    ptm <- proc.time()
    fit = sampling(
      attr(sp,'sm'),
      data = sp,
      iter = 200,
      control = list(adapt_delta = 0.99),
      sample_file = tempdir()
    )
    elapsedTime = (proc.time() - ptm)[3]
    ev = evalSP_lp(sp = sp,
                   fit = fit,
                   X = X_test)
    return(mean((ev$y_hat - y_test) ^ 2))
  }

ghs_mses <-
  function(X_train,
           X_test,
           y_train,
           y_test,
           group_ids,
           num_groups) {
    n = ncol(X_train)
    
    p0s =  c(0.1, 0.25, 0.5) * n
    tau0s = p0s / (sqrt(n) * (n - p0s))
    mses = rep(NA, length(p0s))
    
    i = 1
    for (tau0 in tau0s) {
      mses[i] = ghs_mse(X_train,
                        X_test,
                        y_train,
                        y_test,
                        group_ids,
                        num_groups,
                        tau0)
      i = i + 1
    }
    
    return(mses)
  }

gl_mse <- function(X_train,
                   X_test,
                   y_train,
                   y_test,
                   group_ids) {
  
  cv <- cv.gglasso(
    x = X_train,
    y = y_train,
    group = group_ids,
    loss = "ls",
    pred.loss = "L2",
    lambda.factor = 0.05,
    nfolds = 5
  )
  
  m = gglasso::gglasso(
    x = X_train,
    y = y_train,
    group = group_ids,
    loss = "ls",
    lambda = cv$lambda.min
  )
  return(mean((predict(m, X_test) - y_test) ^ 2))
}

ss_mse <- function(xtrain, train, xtest, y_test) {
  formula = formula(paste("price ~ ", paste(colnames(xtrain), collapse = " + "), sep = ""))
  m <- spikeSlabGAM(formula, data = train)
  
  tmp = sapply(xtest,factor)
  mean((y_test - predict(m, tmp)) ^ 2)
}

lm_mse <- function(train, test, xtest) {
  l = lm(price ~ ., train)
  return(mean((predict(l, xtest) - test$price) ^ 2))
}

computeMSEs <- function(ss, am){
  train = am[ss == 1,]
  test = am[ss == 2,]
  
  xtrain = train %>% select(-"price")
  xtest  = test %>% select(-"price")
  
  mm_train = model.matrix(price ~ . , train)
  
  group_ids = attr(mm_train, "assign")[2:(ncol(mm_train) - 1)]
  num_groups = length(unique(group_ids))
  X_train = mm_train[, 2:(ncol(mm_train) - 1)]
  
  mm_test = model.matrix(price ~ . , test)
  X_test = mm_test[, 2:(ncol(mm_test) - 1)]
  
  y_train = train$price
  y_test  = test$price
  
  return(list(
    ghs_mses = ghs_mses(X_train,X_test, y_train, y_test, group_ids, num_groups),
    gl_mse = gl_mse(X_train, X_test, y_train, y_test, group_ids),
    ss_mse = ss_mse(xtrain, train, xtest, y_test),
    lm_mse = lm_mse(train, test, xtest)
  ))
}

doShuffleSplit <- function(ss) {
  results = unlist(lapply(data.frame(ss$ss_design), function(x)
    computeMSEs(x, ss$am)))
  
  cnames = c("p1","p2","p3","gl","ss","lm")
  
  dfres = data.frame(matrix(
    as.vector(results),
    ncol = length(cnames),
    byrow = T
  ))
  colnames(dfres) = cnames
  
  return(list(results = results,
              sresults = dfres))
}

permute <- function(x){
  sample(x,length(x),replace = F)
}

if (TRUE) {
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  am = downloadAM(
    'https://archive.ics.uci.edu/ml/machine-learning-databases/autos/imports-85.data'
  ) %>% prepareAM()
  
  set.seed(1)    
  ptrain = 0.8
  repeats = 50
  sparse_am = am %>% mutate_if(is.factor,list('p' = ~permute(.)))
  very_sparse_am = sparse_am %>% mutate_if(is.factor,list('p' = ~permute(.)))
  
  ss = createShuffleSplits(am, ptrain = ptrain, repeats = repeats)
  sparse_ss = createShuffleSplits(sparse_am, ptrain = ptrain, repeats = repeats)
  very_sparse_ss = createShuffleSplits(very_sparse_am, ptrain = ptrain, repeats = repeats)
  
  results = doShuffleSplit(ss)
  save(results, file = "automobile.RData")
  sparse_results = doShuffleSplit(sparse_ss)
  save(sparse_results, file = "sparse_automobile.RData")
  very_sparse_results = doShuffleSplit(very_sparse_ss)
  save(very_sparse_results, file = "very_sparse_automobile.RData")  
}

