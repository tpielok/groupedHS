library("splines")

diffMatrix <- function(d, n){
  D = matrix(0, n-d, n)
  for(i in 0:d){
    D = D + diag(n+i)[(i+1):(n-d+i),1:n]*choose(d,i)*(-1)^((d+1)%/%2 + i)
  }
  return(D)
}

penaltyMatrix <- function(d,n){
  #todo vector edge case
  D = diffMatrix(d,n)
  return (t(D) %*% D)
}

dimParams <- function(Blist){
  ret = vector(mode = "numeric", length=length(Blist)+2)
  i = 3
  ret[2] = 1
  for(B in Blist){
       ret[2] = ret[2] * dim(B)[2] 
       ret[i] = dim(B)[2]
       i = i + 1
       ret[1] = dim(B)[1]
  }
  return(ret)
}

BMatrix <- function(Blist){
  dims = dimParams(Blist)
  counter = rep(1,length(dims)-2)
  BM = matrix(1, nrow = dims[1], ncol=dims[2])
  for(i in 1:dims[1]){
    for(j in 1:dims[2]){
      cnum = 1
      for(B in Blist){
        BM[i,j] = BM[i,j] * B[i,counter[cnum]]
        cnum = cnum + 1
      }
      
      toAdd = T
      for(c in 1:length(counter)){
        if(toAdd){
          counter[c] = counter[c] +1
          toAdd = F
        }
        if(counter[c] > dims[2+c]){
          counter[c] = 1
          toAdd =T
        }
      }
    }
  }
  return(BM) 
}

psplineND <- function(Blist, ds){
  dimSpline = length(Blist)
  num_param1d = vector(mode = "numeric", length=dimSpline)
  
  K = vector(mode = "numeric", length=0)
  for(i in 1:dimSpline){
    num_param1d[i] = ncol(Blist[[i]])
    K = c(K,as.vector(penaltyMatrix(ds[i],num_param1d[i])))
  }
  BM <- BMatrix(Blist)
  num_paramnd <- ncol(BM)
  B <- as.vector(BM)
  
  rankK = sum(num_param1d - ds)
  
  pspline <- list(
    dimSpline = dimSpline,
    num_param1d = num_param1d,
    K = K,
    num_paramnd = num_paramnd,
    B = B,
    rankK = rankK
  )
  class(pspline) <- "psplineND"
  return(pspline)
}

stan_params <- function(Y, scale_global, psplines=list(),X=matrix(nrow=1,ncol=0), group_ids=c(0), num_groups=0 ){
  num_data <- length(Y); 
  num_splines = length(psplines)
  if(num_splines > 0){
    hasNoSF = 0
    K = vector("numeric", length = 0)
    B = vector("numeric", length = 0)
    len_K = 0
    len_B = 0
    dimSplines = vector("numeric", length = num_splines)
    num_param1d = vector("numeric", length = 0)
    num_paramnd = vector("numeric", length = num_splines)
    len_p1d = 0
    rankK = vector("numeric", length = num_splines)
    startB = vector("numeric", length = num_splines)
    startG = vector("numeric", length = num_splines)
    startK = vector("numeric", length = num_splines)
    endK = vector("numeric", length = num_splines)
    startP1d = vector("numeric", length = num_splines)
    num_params = 0
    i = 1
    for(p in psplines){
      K = c(K,p$K)
      len_K = len_K + length(p$K)
      B = c(B,p$B)
      len_B = len_B + length(p$B)
      if(i == 1){
        startB[1] = 1
        startG[1] = 1
        startK[1] = 1
        startP1d[1] = 1
        endK[1] = length(p$K)
      }else{
        presp = psplines[[i-1]]
        startB[i] = startB[i-1] + length(presp$B)
        startG[i] = startG[i-1] + presp$num_paramnd
        startK[i] = startK[i-1] + length(presp$K)
        startP1d[i] = startP1d[i-1] + dimSplines[i-1]
      }
      endK[i] = startK[i] + length(p$K)-1
      
      dimSplines[i] = p$dimSpline
      num_param1d = c(num_param1d, p$num_param1d)
      len_p1d = len_p1d + length(p$num_param1d)
      num_paramnd[i] = p$num_paramnd
      num_params = num_params + p$num_paramnd
      rankK[i] = p$rankK
      i = i + 1
    }
  }
  else{
    num_splines = 1
    hasNoSF = 1
    K = vector("numeric", length = 1)
    B = vector("numeric", length = 1)
    len_K = 1
    len_B = 1
    dimSplines = vector("numeric", length = 1)
    num_param1d = vector("numeric", length = 1)
    num_paramnd = vector("numeric", length = 1)
    len_p1d = 1
    rankK = vector("numeric", length = 1)
    startB = vector("numeric", length = 1)
    startG = vector("numeric", length = 1)
    startK = vector("numeric", length = 1)
    endK = vector("numeric", length = 1)
    startP1d = vector("numeric", length = 1)
    num_params = 1
  }
  
  num_linparam = ncol(X)
  if(num_linparam == 0){
    num_linparam = 1
    hasNoLP = 1
    num_groups = 1
    X = matrix(1,num_data,1)
    group_ids = matrix(1,1,1)
  }
  else{
    hasNoLP = 0
    group_ids = matrix(group_ids)
  }
  
  return(list(
    num_data = num_data,
    Y = Y,
    scale_global = scale_global,
    
    hasNoLP = hasNoLP,
    num_linparam = num_linparam,
    num_groups = num_groups,
    group_ids = group_ids,
    X = X,
    
    hasNoSF = hasNoSF,
    num_splines = num_splines,
    len_K = len_K,
    len_B = len_B,
    startB = as.matrix(startB),
    startG = as.matrix(startG),
    startK = as.matrix(startK),
    startP1d = as.matrix(startP1d),
    endK = as.matrix(endK),    
    len_p1d = len_p1d, 
    K = as.matrix(K),
    B = as.matrix(B),
    num_params = num_params,
    dimSplines = as.matrix(dimSplines),
    num_param1d = as.matrix(num_param1d),
    num_paramnd = as.matrix(num_paramnd),
    rankK = as.matrix(rankK)
  ))
}

