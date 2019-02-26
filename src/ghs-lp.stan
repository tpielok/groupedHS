data { 
  int num_data;                // number of observations
  vector[num_data] Y;          // outcome variable
  real<lower=0> scale_global;  // scale for tau
  int<lower=0> num_linparam;   // overall number of predictor coefficients
  int<lower=0> num_groups;     // number of groups
  int<lower=0> group_ids[num_linparam,1]; // assignments of levels 
  matrix[num_data, num_linparam] X; // model matrix of the inputs
  real<lower=1> nu_global;     // degree of freedom for the half-t prior for tau
  real<lower=1> nu_local;      // degree of freedom for the half-t prior for lambdas
} 
 
parameters { 
  real beta0;                           // intercept
  real logsigma;                        // log of noise std
  // auxiliary variables 
  vector[num_linparam] z;               
  real<lower=0> r1_global;
  real<lower=0> r2_global;
  vector<lower=0>[num_groups] r1_localB; 
  vector<lower=0>[num_groups] r2_localB;    
} 
 
transformed parameters { 
  real sigma;                                        // noise std
  real<lower=0> tau;                                 // global shrinkage parameter
  vector<lower=0>[num_groups] lambdaB;               // local shrinkage parameter
  vector[num_linparam] beta;                         // regression coefficients
  vector[num_data] Y_hat = rep_vector(0, num_data);  // eta
  sigma = exp(logsigma);
  tau = r1_global * sqrt(r2_global);
  lambdaB = r1_localB .* sqrt(r2_localB);
  for (m in 1:num_linparam) {
      beta[m] = z[m]* lambdaB[group_ids[m,1]]*tau ;
  }
  Y_hat = Y_hat + X*beta + beta0;
}   
 
model {
  r1_global ~ normal(0.0, scale_global*sigma);
  r2_global ~ inv_gamma(0.5*nu_global, 0.5*nu_global);
  r1_localB ~ normal(0.0, 1);
  r2_localB ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
  
  z ~ normal(0, 1);
  Y ~ normal(Y_hat, sigma); 
}
