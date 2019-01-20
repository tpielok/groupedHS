## Evalute lp results

load("lp-results.RData")

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
    boxplot(matrix(mses, ncol=length(p0)),
            ylab=as.character(res$sp),
            xlab=as.character(res$SNR))
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
    boxplot(matrix(rmses, ncol=2*length(p0)),
            ylab=as.character(res$sp),
            xlab=as.character(res$SNR))
  }
}


for (cnt_sp in 1:length(sp_level)){
  for (cnt_SNR in 1:length(SNR)) {
    rmses = vector("numeric", 0)
    for (cnt_p0 in 1:length(p0)){
      tmp1 = vector("numeric", num_innerLoops)
      tmp2 = vector("numeric", num_innerLoops)
      for (ii in 1:num_innerLoops) {
        res = results[[cnt_sp]][[cnt_SNR]][[cnt_p0]][[ii]]
        tmp1[ii] = res$rmse_levels[1,3]
        tmp2[ii] = res$rmse_levels[2,3]
      }
      rmses = c(rmses, tmp1,tmp2)
    }
    boxplot(matrix(rmses, ncol=2*length(p0)),
            ylab=as.character(res$sp),
            xlab=as.character(res$SNR))
  }
}


