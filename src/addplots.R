library("latex2exp")
library(ggplot2)

dshrink <- function(kappa, tau=1,sigma=1,n=1){
  f = sigma^-1*sqrt(n)*tau
  return(f/pi/((f^2-1)*kappa+1)/sqrt(kappa)/sqrt(1-kappa))
}


x =  c(seq(0, 0.01, length=10000),seq(0.01, 0.99, length= 100),
       seq(0.99,1,length=10000))#
x = seq(0,1,length.out = 10000)
hx1 = dshrink(x)
hx2 = dshrink(x,tau=0.1)

pdf("../plots/kappa.pdf")
plot(x,hx1,type = "l", xlab = TeX('$\\kappa_j$'), ylab = "Density",
     ylim=c(0,3), xlim=c(0,1), cex=3)
lines(x,hx2,lty=29)
df = data.frame(x, hx1, hx2)

ggplot(data = data.frame(x), aes(x)) + 
  ylim(0,2) + 
  xlim(0,1) +
  stat_function(fun = dshrink) + 
  stat_function(fun = function(x) dshrink(x,tau=0.1), linetype=2) 

ggplot(data=df, aes(x,hx1)) + 
  geom_line() + 
  geom_point() + 
  geom_line(aes(x,hx2),linetype=2) + 
  xlim(0,1) + 
  ylim(0,3)

dev.off()

