library("rstan")
source("./ghsAM.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
sm = stan_model("ghsAM.stan")

resfile = "lp-results.RData"

n = 100 # to tune
n_train = 50

glevels = c(3, 5, 9)
num_dg  = length(glevels)

sp_level = c(1, 2) # Level of sparsity

#SNR = c(0.1,0.5, 1, 5,10) # to tune

SNR = c(0.1, 1, 5)

beta_min = 1
beta_max = 2

num_innerLoops = 100

templ_betas = vector("numeric", 0)
templ_group_ids = vector("numeric", 0)

i = 1
for (gl in glevels) {
  templ_betas = c(templ_betas, seq(beta_min, beta_max, length.out = gl - 1))
  templ_group_ids = c(templ_group_ids, rep(i, gl - 1))
  i = i + 1
}

num_p0 = 3

#results = l <- vector("list",num_p0 * length(sp_level) * length(SNR) * num_innerLoops)
seed_cnt = 1

results = vector("list", length(sp_level))


for (cnt_sp in 1:length(sp_level)) {
  cur_sp_level = sp_level[cnt_sp]
  
  num_nonzero = length(templ_betas)
  betas = vector("numeric", 0)
  group_ids = vector("numeric", 0)
  
  betas = c(templ_betas, rep(0, num_nonzero * cur_sp_level))
  
  p0 = c(1, num_nonzero, length(betas) - 1)
  
  group_ids = templ_group_ids
  for (i in 1:cur_sp_level) {
    group_ids = c(group_ids, templ_group_ids + i * num_dg)
  }
  
  num_groups = (cur_sp_level + 1) * num_dg
  
  results[[cnt_sp]] = vector("list", length(SNR))
  
  for (cnt_SNR in 1:length(SNR)) {
    cur_SNR = SNR[cnt_SNR]
    results[[cnt_sp]][[cnt_SNR]] = vector("list", length(p0))
    for (cnt_p0 in 1:length(p0)){
      cur_p0 = p0[cnt_p0]
      results[[cnt_sp]][[cnt_SNR]][[cnt_p0]] = vector("list", num_innerLoops)
      for (ii in 1:num_innerLoops) {
        
        set.seed(seed_cnt)
        inpMat = matrix(nrow = n, ncol = num_dg * (1 + cur_sp_level))
        for (i in 1:ncol(inpMat)) {
          inpMat[, i] = as.character(sample(1:glevels[((i - 1) %% num_dg) + 1], n, TRUE))
        }
        
        df = data.frame(inpMat)
        mM = model.matrix( ~ ., df)
        X = mM[, 2:ncol(mM)] # Intercept gets estimated seperately
        
        eta = X[1:n_train, ] %*% betas
        
        sd_eta = sqrt(n_train / (n_train - 1)) * sd(eta)
        
        Y = eta + rnorm(n_train, 0, sd_eta ^ 2 / cur_SNR)
        
        #df$Y = Y
        
        #lm(Y ~ ., data=df)
        
        tau0 = cur_p0 / (sqrt(n) * (length(betas) - cur_p0))
        s_p = stan_params(
          Y = Y[, 1],
          scale_global = tau0,
          X = X[1:n_train, ],
          group_ids = group_ids,
          num_groups = num_groups
        )
        
        ptm <- proc.time()
        fit = sampling(
          sm,
          data = s_p,
          iter = 500,
          control = list(adapt_delta = 0.99),
          seed = seed_cnt,
          sample_file=tempdir()
        )
        elapsedTime = (proc.time() - ptm)[3]
        
        la = extract(fit, permuted = TRUE)
        
        beta_est = colMeans(la$beta)
        
        res = betas - beta_est
        sq_res = (res) ^ 2
        
        rmse_nonzero = sqrt(mean(sq_res[1:num_nonzero]))
        rmse_zero    = sqrt(mean(sq_res[(num_nonzero + 1):length(sq_res)]))
        
        rmse_levels = matrix(nrow = 2, ncol = num_dg) # 1: nonzero 2: zero
        me_levels =  matrix(nrow = 2, ncol = num_dg) # 1: nonzero 2: zero
        
        sq_res = matrix(sq_res, nrow = cur_sp_level + 1, byrow = TRUE)
        res = matrix(res, nrow = cur_sp_level + 1, byrow = TRUE)
        
        offset = 1
        for (l in 1:num_dg) {
          rmse_levels[1, l] = sqrt(mean(sq_res[1, offset:(offset + glevels[l] - 2)]))
          me_levels[1, l] = mean(res[1, offset:(offset + glevels[l] - 2)])
          
          rmse_levels[2, l] =  sqrt(mean(colMeans(
            matrix(sq_res[2:(cur_sp_level + 1), offset:(offset + glevels[l] - 2)], nrow =
                     cur_sp_level)
          )))
          me_levels[2, l] = mean(colMeans(matrix(res[2:(cur_sp_level + 1), offset:(offset +
                                                                                     glevels[l] - 2)], nrow = cur_sp_level)))
          offset = offset + glevels[l] - 1
        }
        
        rmse_levels
        me_levels
        
        Y_test = X[n_train:n, ] %*% betas
        Y_est = X[n_train:n, ] %*% beta_est
        
        beta0_est = mean(la$beta0)
        sigma_est = mean(la$sigma)
        
        mean_y = mean(Y)
        sd_y   = sd(Y)
        mse    = sqrt(mean((Y_test - Y_est) ^ 2))
        result = list(
          "rmse_nz" = rmse_nonzero,
          "rmse_z" = rmse_zero,
          "rmse_levels" = rmse_levels,
          "me_levels" = me_levels,
          "mse" = mse,
          "mean_y" = mean_y,
          "sd_y" = sd_y,
          "beta0_est" = beta0_est,
          "beta_est" = beta_est,
          "sigma_est" = sigma_est,
          "sp" = cur_sp_level,
          "SNR" = cur_SNR,
          "p0" = cur_p0,
          "tau0" = tau0,
          "seed" = seed_cnt,
          "elapsedTime" = elapsedTime,
          "ii" = ii,
          "seed_cnt" = seed_cnt
        )
        
        results[[cnt_sp]][[cnt_SNR]][[cnt_p0]][[ii]] = result
        print(seed_cnt/(length(p0)*length(SNR)*length(sp_level)*num_innerLoops))
        seed_cnt = seed_cnt + 1
      }
      save(results, file = resfile)
    }
  }
}

