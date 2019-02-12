library("rstan")
library("rgl")
source("./ghsAM.R")
source("./Util.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(123)
X <- matrix(runif(n = 300, min = -5, max = 5),ncol=3)

Y_true <- sin(X[,1]) + cos(X[,2])

plot3d(x=X[,1],y=X[,3],z=Y_true,
       col = map2color(Y_true, heat.colors(20)))

Y <- Y_true + rnorm(length(Y_true),0,.2) # adding noise

B1 <- bs(X[,1], knots=seq(-5,5,2), degree=2, intercept = TRUE) # creating the B-splines
B2 <- bs(X[,2], knots=seq(-5,5,2), degree=3, intercept = TRUE)

d1 = c(3)
p1 <- psplineND(list(B1),d1)

d2 = c(3)
p2 <- psplineND(list(B2),d2)


s_p <- stan_params(Y,list(p1,p2), scale_global = 1)

s_p

sm<-stan_model("ghsAM.stan") 
fit<-sampling(sm, data = s_p,iter=500,control = list(adapt_delta=0.99999))

la <- extract(fit, permuted = TRUE)

length(X)
mean(la$gamma)
length(la$Y_hat)

gamma = colMeans(la$gamma)
beta0 = mean(la$beta0)

length(gamma)

B = matrix(p1$B, nrow=length(Y))
y_hat = B  %*% gamma[1:ncol(B)]
offset = ncol(B)
plot(X[,1], y_hat)

B = matrix(p2$B, nrow=length(Y))
y_hat = B  %*% gamma[(offset+1):(offset+ncol(B))]

plot(X[,2], y_hat)

y_hat = beta0 + matrix(p1$B, nrow=length(Y)) %*% gamma
