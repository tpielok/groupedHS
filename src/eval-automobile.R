library("ggplot2")

load(file = "data/automobile.RData")
load(file = "data/sparse_automobile.RData")
load(file = "data/very_sparse_automobile.RData")

nc = ncol(results$sresults)
nr = nrow(results$sresults)
l = nc*nr
df = data.frame(matrix(NA,nrow = 3*l, ncol=3))
colnames(df) = c("mse","method","sparse")
df$mse[1:l] = as.vector(as.matrix(sqrt(results$sresults)))
df$mse[(l+1):(2*l)] = as.vector(as.matrix(sqrt(sparse_results$sresults)))
df$mse[(2*l+1):(3*l)] = as.vector(as.matrix(sqrt(very_sparse_results$sresults)))
df$method = rep(as.vector(sapply(as.character(1:nc), function(x) rep(x,nr))),3)
df$sparse = c(rep("1",l),rep("2",l),rep("3",l))

df$method = factor(df$method, labels=c("p0 = 0.1D","p0 = 0.25D","p0 = 0.5D", "gglasso","spikeSlabGAM","linear model"))
df$sparse = factor(df$sparse, labels=c("unaltered","heightened (1x)","heightened (2x)"))

ggplot(data = df, aes(y=mse,x=method)) + 
  geom_boxplot(aes(fill = sparse, color=sparse)) + 
  xlab("method") + ylab("RMSE") +
  theme(legend.position = "bottom",axis.text.x=element_text(angle=45, hjust=1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12)) +
  scale_fill_manual(name = "Level of sparsity",
                    values = c("grey40", "lightskyblue","chartreuse2")) +
  scale_color_manual(name ="Level of sparsity", values = c("black", "blue","forestgreen")) 

ggsave("../plots/am-mse.pdf")


