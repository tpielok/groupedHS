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
} 
 
parameters { 
  vector[num_params] gamma; 
  real beta0; 
  real<lower=0> sigma; 

  vector<lower=0>[num_splines] lambdaG; 
  vector<lower=0>[num_groups] lambdaB; 
  vector[num_linparam-hasNoLP] z;
  real<lower=0> tau; 
} 
 
transformed parameters { 
  vector[num_linparam-hasNoLP] beta;
  vector[num_data] Y_hat = rep_vector(0, num_data); 
  
  if(hasNoSF == 0){
    for (i in 1:num_splines){
      Y_hat = 
        Y_hat + to_matrix(B[
          startB[i,1]:(startB[i,1]+num_data*num_paramnd[i,1]-1),1],
          num_data, num_paramnd[i,1]) *
          gamma[startG[i,1]:(startG[i,1]+num_paramnd[i,1]-1)];
    }
  }
  
  if(hasNoLP == 0){
    for (m in 1:num_linparam) {
      beta[m] = z[m]* lambdaB[group_ids[m,1]]*tau ;
  	 }
    Y_hat = Y_hat + X*beta;
  }
  
  Y_hat = Y_hat + beta0;
}   
 
model {
  sigma ~ normal(0, 1);  
  tau ~ cauchy(0, scale_global*sigma);
  
  if(hasNoSF == 0){  
    lambdaG ~ cauchy(0, 1);
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
    z ~ normal(0, 1);
    lambdaB ~ cauchy(0, 1);
  }

  Y ~ normal(Y_hat, sigma); 
}
