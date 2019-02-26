data { 
  int<lower=0> num_data;              // number of observations
  vector[num_data] Y;                 // outcome
  int<lower=0> num_splines;           // number of splines
  int<lower=0> num_params;            // number of predictor coefficients
  int<lower=0> len_Bs;                // length of flattened transformed model model matrices Bs
  int<lower=0> bStart[1, num_splines];// starting indices of singular B
  int<lower=0> bEnd[1, num_splines];  // ending indices of singular B
  vector[len_Bs] Bs;                  // flattened transformed model model matrices Bs
  int<lower=0> gStart[1, num_splines];// starting indices of coefficients of a singluar predictor 
  int<lower=0> gLen[1, num_splines];  // number of coefficients of a singluar predictor 
  int<lower=0> gEnd[1, num_splines];  // ending indices of coefficients of a singluar predictor 
  
  real<lower=1> nu_global;            // degree of freedom for the half-t prior for tau
  real<lower=1> nu_local;             // degree of freedom for the half-t prior for u
  real<lower=0> scale_global;         // scale for tau
} 
 
parameters { 
  real beta0;                         // intercept
  real logsigma;                      // log noise std
  // auxiliary variables
  vector[num_params] z; 
  real<lower=0> r1_global;
  real<lower=0> r2_global;
  vector<lower=0>[num_splines] r1_localG; 
  vector<lower=0>[num_splines] r2_localG; 
} 
 
transformed parameters { 
  real sigma;                        // noise std
  real tau;                          // global shrinkage parameter
  vector[num_splines] lambda;        // local shrinkage parameter
  
  vector[num_params] u;              // coefficients of all predictors
  vector[num_data] Y_hat = rep_vector(0, num_data); // eta
  tau = r1_global * sqrt(r2_global);
  lambda = r1_localG .* sqrt(r2_localG);
  
  sigma = exp(logsigma);
  
  for (i in 1:num_splines){
      u[gStart[1,i]:gEnd[1,i]] = z[gStart[1,i]:gEnd[1,i]] * lambda[i]* tau;
      Y_hat = Y_hat + to_matrix(Bs[bStart[1,i]:bEnd[1,i]],
          num_data, gLen[1,i]) * u[gStart[1,i]:gEnd[1,i]];
          
  }
  Y_hat = Y_hat + beta0;
}   
 
model {
  r1_global ~ normal(0.0, scale_global*sigma);
  r2_global ~ inv_gamma(0.5*nu_global, 0.5*nu_global);
  
  r1_localG ~ normal(0.0, 1);
  r2_localG ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
  
  z ~ normal(0,1);
  Y ~ normal(Y_hat, sigma); 
}