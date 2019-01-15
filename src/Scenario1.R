library("rstan")
library("dplyr")
source("./ghsAM.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(123)

# Data creation
num_data = 50
num_cont = 15
num_dscrt = 15

max_lev = 6

num_groups = num_cont + num_dscrt
group_ids = 1:num_cont

num_lev = rep(1, num_groups)
inpMat = matrix(nrow = num_data, ncol=num_groups)

inpMat[1:num_data, 1:num_cont] = matrix(runif(num_data*num_cont), nrow = num_data, ncol = num_cont)
for(i in 1:num_dscrt){
  num_lev[num_cont+i] = sample(2:max_lev,1)
  group_ids = c(group_ids, rep(i,num_lev[num_cont+i]-1))
  inpMat[1:num_data, num_cont+i] = sample(1:num_lev[num_cont+i], num_data, TRUE)
}

df = data.frame(inpMat)
for(i in num_groups:(num_cont+1)){
  df[,i] = factor(df[,i])
}

mM = model.matrix(~ ., df)

X = mM[,2:ncol(mM)] # Intercept gets estimated seperately
group_ids
beta_real= c(5, -2, 1, -3)
beta_0 = -2
Y = beta_real[1]*X[,1] + beta_real[2]*X[,2] + beta_real[3]*X[,16] + beta_real[4]*X[,17] +
  rnorm(num_data, 0, 0.2) + beta_0

df$Y = Y

lm(Y ~ ., data=df)

p0 = length(beta_real) - 1
tau0 = p0/(sqrt(num_data)*(ncol(X) - p0))

s_p = stan_params(Y = Y,scale_global = tau0, X=X, group_ids = group_ids, num_groups = num_groups)

sm = stan_model("ghsAM.stan") 
fit = sampling(sm, data = s_p,iter=500,control = list(adapt_delta=0.99))

la = extract(fit, permuted = TRUE)

pos_betas = as_tibble(la$beta)

hs_betas = pos_betas %>% 
  select_if(function(col) 
    sign(quantile(col,probs = 0.05)) ==
      sign(quantile(col,probs = 0.95))
  ) %>% 
  mutate(V0 = la$beta0) %>% 
  summarise_all(funs(mean)) 

hs_betas


