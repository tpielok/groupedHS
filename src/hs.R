library("dplyr")
library("rstan")
library("glmnet")
# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Create data
set.seed(1234)
n <- 50
d <- 60

beta_real = c(3, -3, 2, -1)
var_real = 0.5

inp_raw <- runif(n*d)
inp_mat <- matrix(inp_raw,
                  nrow = n,
                  ncol = d)

data <- as_tibble(inp_mat)
data <- mutate(data,
               y = beta_real[1] + beta_real[2]*V3 +  beta_real[3]*V4  +  beta_real[4]*V5 + rnorm(n, 0, var_real))

# Horseshoe

p0 <- length(beta_real) - 1
tau0 <- p0/(sqrt(n)*(d - p0))

stan_params <- list(
  n = n,
  d = d,
  y = data$y,
  x = inp_mat,
  
  scale_icept = 100,
  scale_global = tau0
)


fit <- stan(file = 'hs.stan', data = stan_params, 
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

# Lasso

x <- model.matrix(y~., data)[,-1]
y <- data$y
lambda <- 10^seq(10, -2, length = 100)

cv.out <- cv.glmnet(x, y, alpha = 1)
bestlam <- cv.out$lambda.min
lasso.mod <- glmnet(x, y, alpha = 1, lambda = bestlam)
lasso_betas <- as_tibble(as.matrix(t(coef(lasso.mod)))) %>% select_if(. != 0)

# Estimated selected variables

print(beta_real)
print(hs_betas)
print(lasso_betas)

