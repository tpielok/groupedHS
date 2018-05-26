library("dplyr")
library("rstan")
#library("glmnet")
library("gglasso")
# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Create data
set.seed(1234)
n <- 50
dngrp <- 10
dgrp <- 6
dlev <- 3
d <- dngrp + dgrp * dlev

num_groups <- dngrp + dgrp
group <- c(1:dngrp,rep((dngrp+1):(dngrp+dgrp),each=dlev))

beta_real = c(3, -3, 2, -1, 3, -2.5)
var_real = 0.5

inp_raw <- runif(n*dngrp)
inp_raw <- c(inp_raw, rbinom(dgrp*dlev*n, 1, 0.5))
inp_raw <- numeric()
for (i in seq(1:n)){
  inp_raw <- c(inp_raw, runif(dngrp))
  for (j in seq(1:dgrp)){
    tmp <- rep(0, dlev)
    tmp[sample(1:dlev,1)] <- 1
    inp_raw <- c(inp_raw, tmp)
  }
}

inp_mat <- matrix(inp_raw,
                  ncol = d,
                  byrow = T)

data <- as_tibble(inp_mat)
data <- mutate(data,
               y = beta_real[1] + beta_real[2]*V3 +  
                 beta_real[3]*V4  +  beta_real[4]*V11 +
                 beta_real[5]*V12 + beta_real[6]*V13 +
                 rnorm(n, 0, var_real))
# Grouped Horseshoe
#tau0 <- 1
p0 <- length(beta_real) - 1
tau0 <- p0/(sqrt(n)*(d - p0)) # todo


stan_params <- list(
  n = n,
  d = d,
  y = data$y,
  x = inp_mat,
  g = num_groups,
  group_ids = group,
  scale_icept = 100,
  scale_global = tau0
)

fit <- stan(file = 'ghs.stan', data = stan_params, 
            iter = 1000, chains = 4, control = list(adapt_delta = 0.99))

la <- extract(fit, permuted = TRUE)
pos_betas <- as_tibble(la$beta)

hs_betas <- pos_betas %>% 
  select_if(function(col) 
    sign(quantile(col,probs = 0.05)) ==
      sign(quantile(col,probs = 0.95))
  ) %>% 
  mutate(V0 = la$beta0) %>% 
  summarise_all(funs(mean)) 

# Grouped Lasso



# Results
print(hs_betas)
