## Evalute lp results

load("lp-results.RData")

p0 = c(1,2,3)

par(mfrow=c(length(sp_level),length(SNR)))

for (cnt_sp in 1:length(sp_level)){
  for (cnt_SNR in 1:length(SNR)) {
    mses = vector("numeric", 0)
    for (cnt_p0 in 1:length(p0)){
      tmp = vector("numeric", num_innerLoops)
      for (ii in 1:num_innerLoops) {
        res = results[[cnt_sp]][[cnt_SNR]][[cnt_p0]][[ii]]
        tmp[ii] = res$mse
      }
      mses = c(mses, tmp)
    }
    
    if(res$sp == 1){
      yl = "lo"
    }else{
      yl = "high"
    }
   
    boxplot(matrix(mses, ncol=length(p0)),
            ylab=yl,
            xlab=paste("SNR",as.character(res$SNR)),
            names=c("p1","p2","p3"))
    mtext("MSE")
  }
}



par(mfrow=c(length(sp_level),length(SNR)))

for (cnt_sp in 1:length(sp_level)){
  for (cnt_SNR in 1:length(SNR)) {
    rmses = vector("numeric", 0)
    for (cnt_p0 in 1:length(p0)){
      tmp1 = vector("numeric", num_innerLoops)
      tmp2 = vector("numeric", num_innerLoops)
      for (ii in 1:num_innerLoops) {
        res = results[[cnt_sp]][[cnt_SNR]][[cnt_p0]][[ii]]
        tmp1[ii] = res$rmse_nz
        tmp2[ii] = res$rmse_z
      }
      rmses = c(rmses, tmp1,tmp2)
    }

    if(res$sp == 1){
      yl = "lo"
    }else{
      yl = "high"
    }
    
    boxplot(matrix(rmses, ncol=2*length(p0)),
            ylab=yl,
            xlab=paste("SNR",as.character(res$SNR)),
            names=c("p1NZ","p1Z","p2NZ","p2Z","p3NZ","p1Z"))
    mtext("Overall quality (param mse)")    
  }
}



for(i in 1:3){

for (cnt_sp in 1:length(sp_level)){
  for (cnt_SNR in 1:length(SNR)) {
    rmses = vector("numeric", 0)
    for (cnt_p0 in 1:length(p0)){
      tmp1 = vector("numeric", num_innerLoops)
      tmp2 = vector("numeric", num_innerLoops)
      for (ii in 1:num_innerLoops) {
        res = results[[cnt_sp]][[cnt_SNR]][[cnt_p0]][[ii]]
        tmp1[ii] = res$rmse_levels[1,i]
        tmp2[ii] = res$rmse_levels[2,i]
      }
      rmses = c(rmses, tmp1,tmp2)
    }
    if(res$sp == 1){
      yl = "lo"
    }else{
      yl = "high"
    }
    
    if(i==1){
      si = "small"
    }else if(i==2){
      si = "medium"
    }else{
      si = "big"
    }
    
    boxplot(matrix(rmses, ncol=2*length(p0)),
            ylab=yl,
            xlab=paste("SNR",as.character(res$SNR)),
            names=c("p1NZ","p1Z","p2NZ","p2Z","p3NZ","p1Z"))
    mtext(paste(si,"group (param mse)"))
  }
}

}
