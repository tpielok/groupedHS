#include /psplineFunc.stan

data { 
  int num_data;
  vector[num_data] Y;   
  real<lower=0> scale_global;
    
  // Linear terms
  int<lower=0> hasNoLP;
  int<lower=0> num_linparam; 
  int<lower=0> num_groups; 
     
  int<lower=0> group_ids[num_linparam,1];
  matrix[num_data, num_linparam] X;
  
  // Smooth Funcions
  //int<lower=0> hasNoSF;  
  int<lower=0>  num_splines;
  int<lower=0> hasNoSF;
  int<lower=0> len_K;
  int<lower=0> len_B;
  int<lower=0> len_p1d;
  matrix[len_K,1] K;
  matrix[len_B,1] B;
  int<lower=0> num_params;
  int<lower=0> endK[num_splines,1];
  int<lower=0> startG[num_splines,1];
  int<lower=0> startB[num_splines,1];
  int<lower=0> startK[num_splines,1];      
  int<lower=0> startP1d[num_splines,1];
  int<lower=0> dimSplines[num_splines,1];
  int<lower=0> num_param1d[len_p1d,1];
  int<lower=0> num_paramnd[num_splines,1];
  int<lower=0> rankK[num_splines,1];
  
  real<lower=1> nu_global;
  real<lower=1> nu_local;
  
  //matrix[num_groups,1] scale_localB;
} 
 
parameters { 
  vector[num_params] gamma; 
  real beta0; 
  real logsigma; 
  
  vector[num_linparam-hasNoLP] z;
  
  real<lower=0> r1_global;
  real<lower=0> r2_global;
  
  vector<lower=0>[num_splines-hasNoSF] r1_localG; 
  vector<lower=0>[num_splines-hasNoSF] r2_localG;   
  vector<lower=0>[num_groups-hasNoLP] r1_localB; 
  vector<lower=0>[num_groups-hasNoLP] r2_localB;    
} 
 
transformed parameters { 
  real sigma;
  vector[num_linparam-hasNoLP] beta;
  vector[num_data] Y_hat = rep_vector(0, num_data); 
  
  real<lower=0> tau; 
  vector<lower=0>[num_splines-hasNoSF] lambdaG; 
  vector<lower=0>[num_groups-hasNoLP] lambdaB; 
  
  sigma = exp(logsigma);
  tau = r1_global * sqrt(r2_global);
  
  if(hasNoSF == 0){
    lambdaG = r1_localG .* sqrt(r2_localG);
    for (i in 1:num_splines){
      Y_hat = 
        Y_hat + to_matrix(B[
          startB[i,1]:(startB[i,1]+num_data*num_paramnd[i,1]-1),1],
          num_data, num_paramnd[i,1]) *
          gamma[startG[i,1]:(startG[i,1]+num_paramnd[i,1]-1)];
    }
  }
  
  if(hasNoLP == 0){
    lambdaB = r1_localB .* sqrt(r2_localB);
    for (m in 1:num_linparam) {
      beta[m] = z[m]* lambdaB[group_ids[m,1]]*tau ;
  	 }
    Y_hat = Y_hat + X*beta;
  }
  
  Y_hat = Y_hat + beta0;
}   
 
model {
  r1_global ~ normal(0.0, scale_global*sigma);
  r2_global ~ inv_gamma(0.5*nu_global, 0.5*nu_global);
  
  if(hasNoSF == 0){  
    r1_localG ~ normal(0.0, 1);
    r2_localG ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
    
    for (i in 1:num_splines){
        gamma[startG[i,1]:(startG[i,1]+num_paramnd[i,1]-1)] ~ 
          imprNormal(
            K[startK[i,1]:endK[i,1],1], 
            tau*lambdaG[i], 
            rankK[i,1], 
            dimSplines[i,1], 
            num_param1d[startP1d[i,1]:(startP1d[i,1]+dimSplines[i,1]-1),1]);    
    }
  }
  if(hasNoLP == 0){
    r1_localB ~ normal(0.0, 1);
    r2_localB ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
    z ~ normal(0, 1);
  }

  Y ~ normal(Y_hat, sigma); 
}
