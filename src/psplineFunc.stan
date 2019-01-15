functions {
  int numSteps(int dimSpline, int[] num_param1d, int i){
    int ret = 1;
    for (j in 1:dimSpline){
      if(j != i){
        ret = ret * num_param1d[j];
      }
    }
    return ret;
  }
  
  int getPos(int dimSpline, int[] num_param1d, int[] pos){
    int cpos = 0;
    int offset = 1;
    for (i in 1:dimSpline){
      cpos = cpos + pos[i]*offset;
      offset = offset * num_param1d[i];
    }
    return cpos + 1;
  }
  
  real K1Dir(vector y, vector K, int n, int offset, int len,
    int dimSpline, int[] num_param1d, int i){
    real ret = 0.0; 
    int toAdd;
    matrix[n, n] Kmat = to_matrix(K[offset:(offset+len-1)], n, n); 
    vector[n] param;

    int pos[dimSpline];
    for(j in 1:dimSpline){
      pos[j] = 0;
    }

    for (j in 1:numSteps(dimSpline, num_param1d, i)){
      for(k in 1:num_param1d[i]){
        pos[i] = k-1;
        param[k] = y[getPos(dimSpline, num_param1d, pos)];
      }
      
      ret = ret + param'*Kmat*param;

      // Update Pos
      toAdd = 1;
      for (l in 1:dimSpline){
        if (toAdd && l != i){
          toAdd = 0;
          pos[l] = pos[l] + 1;
          if (pos[l] >= num_param1d[l]){
            pos[l] = 0;
            toAdd = 1;
          }
        }
      }
    }
    
    return ret;
  }
  
  real imprNormal_lpdf(vector y, vector K, real tau, int rankK, int dimSpline, int[] num_param1d) {
    int n;
    int len;
    int offset = 1;         
    real dsq = 0.0;
    
    for (i in 1:dimSpline) {
      n = num_param1d[i];
      len = n*n;
      dsq = dsq + K1Dir(y, K, n, offset, len, dimSpline, num_param1d, i);
      offset = offset + len;
  	}
    return -log(tau)*rankK -1/(2*tau^2)*dsq;
  }
  
  real imprNormalM_lpdf(vector y, matrix K, real tau, int rankK) {
    return -log(tau)*rankK -1/(2*tau^2)*y'*K*y;
  }
} 
