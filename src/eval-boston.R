library("ggplot2")
library("latex2exp")

load(file = "data/Boston.RData")
load(file = "data/sparseBoston.RData")

nc = ncol(results$sresults)
nr = nrow(results$sresults)
l = nc*nr
df = data.frame(matrix(NA,nrow = 2*l, ncol=3))
colnames(df) = c("mse","method","sparse")
df$mse[1:l] = as.vector(as.matrix(sqrt(results$sresults)))
df$mse[(l+1):(2*l)] = as.vector(as.matrix(sqrt(sparse_results$sresults)))

df$method = rep(as.vector(sapply(as.character(1:nc), function(x) rep(x,nr))),2)
df$sparse = c(rep("1",l),rep("2",l))

df$method = factor(df$method, labels=c("p0 = 0.1D","p0 = 0.2D","p0 = 0.4D", "GAMSEL","spikeSlabGAM","linear model","GAM"))
df$sparse = factor(df$sparse, labels=c("unaltered","heightened"))

ggplot(data = df, aes(y=mse,x=method)) + 
  geom_boxplot(aes(fill = sparse, color=sparse)) + 
  xlab("method") + ylab(TeX('$\\log_{10}$(RMSE)')) +
  scale_y_continuous(breaks = c(10,100)) + 
  coord_trans(y = "log10") + 
  theme(legend.position = "bottom",axis.text.x=element_text(angle=45, hjust=1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12)) +
  scale_fill_manual(name = "Level of sparsity",
                    values = c("grey40", "lightskyblue","chartreuse2")) +
  scale_color_manual(name ="Level of sparsity", values = c("black", "blue","forestgreen")) 

ggsave("../plots/bn-mse.pdf", width = 6)

