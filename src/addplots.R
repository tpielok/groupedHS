library("latex2exp")

dshrink <- function(kappa, tau=1,sigma=1,n=1){
  f = sigma^-1*sqrt(n)*tau
  return(f/pi/((f^2-1)*kappa+1)/sqrt(kappa)/sqrt(1-kappa))
}


x =  c(seq(0, 0.05, length=500),seq(0.05, 0.9, length= 100),
       seq(0.9,1,length=100))
hx1 = dshrink(x)
hx2 = dshrink(x,tau=0.1)

pdf("../plots/kappa.pdf")
plot(x,hx1,type = "l", xlab = TeX('$\\kappa_j$'), ylab = "Density",
     ylim=c(0,3), xlim=c(0,1), cex=3)
lines(x,hx2,lty=29)
dev.off()

