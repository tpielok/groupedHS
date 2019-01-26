functions {
  
  real imprNormalM_lpdf(vector y, matrix P, real tau, real rankK) {
    return -log(tau)*rankK -1/(2*tau^2)*y'*P*y;
  }
  
} 
