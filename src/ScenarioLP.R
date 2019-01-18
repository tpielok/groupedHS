library("rstan")
source("./ghsAM.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(1)

n = 50 # to tune
glevels = c(3,5,9)
num_dg  = length(glevels)

sp_level = c(1,2) # Level of sparsity

SNR = c(5, 10, 20) # to tune

beta_min = 1
beta_max = 2

cur_sp_level = sp_level[2]
cur_SNR = SNR[1]

betas = vector("numeric",0)
group_ids = vector("numeric",0)

i = 1
for (gl in glevels){
  betas = c(betas,seq(beta_min,beta_max,length.out = gl-1))
  group_ids = c(group_ids, rep(i, gl-1))
  i = i +1
}
num_nonzero = length(betas)
betas = c(betas, rep(0,num_nonzero*cur_sp_level))

group_ids_tmp = group_ids
for (i in 1:cur_sp_level){
  group_ids = c(group_ids, group_ids_tmp + i*num_dg)
}

num_groups = (cur_sp_level+1)*num_dg

inpMat = matrix(nrow = n, ncol = num_dg*(1+cur_sp_level))
for (i in 1:ncol(inpMat)){
  inpMat[,i] = as.character(sample(1:glevels[((i-1)%% num_dg) + 1], n, TRUE))
}

df = data.frame(inpMat)
mM = model.matrix(~ ., df)
X = mM[,2:ncol(mM)] # Intercept gets estimated seperately

eta = X %*% betas

sd_eta = sqrt(n/(n-1))*sd(eta)
Y = eta + rnorm(n,0,sd_eta^2/cur_SNR)

df$Y = Y

lm(Y ~ ., data=df)

p0 = num_nonzero
tau0 = p0/(sqrt(n)*(length(betas) - p0))

s_p = stan_params(Y = Y[,1],scale_global = tau0, X=X, group_ids = group_ids, num_groups = num_groups)

sm = stan_model("ghsAM.stan") 
fit = sampling(sm, data = s_p,iter=500,control = list(adapt_delta=0.99))

la = extract(fit, permuted = TRUE)

beta_est = colMeans(la$beta)

sq_res = (betas - beta_est)^2
rmse_nonzero = sqrt(mean(sq_res[1:num_nonzero]))
rmse_zero    = sqrt(mean(sq_res[(num_nonzero+1):length(sq_res)]))

rmse_levels = matrix(nrow = 2, ncol= num_dg) # 1: nonzero 2: zero

sq_res = matrix(sq_res, nrow = cur_sp_level+1, byrow = TRUE)

offset = 1
for (l in 1:num_dg){
  print(offset)
  print((offset+glevels[l]-2))
  rmse_levels[1,l] = sqrt(mean(sq_res[1,offset:(offset+glevels[l]-2)]))
  rmse_levels[2,l] =  sqrt(mean(colMeans(matrix(sq_res[2:(cur_sp_level+1),offset:(offset+glevels[l]-2)],nrow=cur_sp_level))))
  offset = offset + glevels[l] - 1
}

rmse_levels

