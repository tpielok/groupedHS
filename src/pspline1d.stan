#include /pspline1DFunc.stan

data { 
  int<lower=0> num_data;
  vector[num_data] Y;
  int<lower=0> num_splines;
  int<lower=0> num_params;
  int<lower=0> len_Bs;
  int<lower=0> len_Ps;
  int<lower=0> bStart[1, num_splines];
  int<lower=0> bEnd[1, num_splines];
  vector[len_Bs] Bs;
  int<lower=0> pStart[1, num_splines];
  int<lower=0> pEnd[1, num_splines];
  vector[len_Ps] Ps;     
  int<lower=0> gStart[1, num_splines];
  int<lower=0> gLen[1, num_splines];  
  int<lower=0> gEnd[1, num_splines];  
  
  int<lower=0> pRanks[1, num_splines];  
  
  real<lower=1> nu_global;
  real<lower=1> nu_local;
  real<lower=0> scale_global;
} 
 
parameters { 
  vector[num_params] gamma; 
  
  real<lower=0> r1_global;
  real<lower=0> r2_global;
  
  vector<lower=0>[num_splines] r1_local; 
  vector<lower=0>[num_splines] r2_local; 
  
  real beta0; 
  real logsigma;
} 
 
transformed parameters { 
  real sigma;
  real tau;
  
  vector[num_splines] lambda;   
  vector[num_data] Y_hat = rep_vector(0, num_data); 

  tau = r1_global * sqrt(r2_global);
  lambda = r1_local .* sqrt(r2_local);
  
  sigma = exp(logsigma);
  
  for (i in 1:num_splines){
      Y_hat = Y_hat + to_matrix(Bs[bStart[1,i]:bEnd[1,i]],
          num_data, gLen[1,i]) *
          gamma[gStart[1,i]:gEnd[1,i]];
  }
  Y_hat = Y_hat + beta0;
}   
 
model {
  r1_global ~ normal(0.0, scale_global*sigma);
  r2_global ~ inv_gamma(0.5*nu_global, 0.5*nu_global);
  
  r1_local ~ normal(0.0, 1);
  r2_local ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
  
  for (i in 1:num_splines){
      gamma[gStart[1,i]:gEnd[1,i]] ~ imprNormalM(
        to_matrix(Ps[pStart[1,i]:pEnd[1,i]],
          gLen[1,i], gLen[1,i]), tau*lambda[i], pRanks[1,i]);  
  }
  
  Y ~ normal(Y_hat, sigma); 
}
