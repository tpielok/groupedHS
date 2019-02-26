library("latex2exp")
library("dplyr")
library("tidyr")
source("util-sf.R")

load(get_sf_file())

scenarios = get_sf_scenarios()

for (i in 1:length(results)){
  r = results[[i]]
  num_cov <- ifelse(scenarios$sparse[r$dgp] == "high", 20, 16)
  num_nz <- ifelse(scenarios$sparse[r$dgp] == "high", 4, 12)
  results[[i]]$qnz = mean(r$quality[1:num_nz])
  results[[i]]$qz = mean(r$quality[(num_nz+1):num_cov])
}


aresults = list()
for (i in 1:nrow(scenarios)){
  aresults[[i]] = list()
}

for (r in results){
  aresults[[r$dgp]][[r$j]] = r
}

# MSE

nsz = nrow(scenarios)
mmses =  matrix(0, nrow = 3 * nsz * 100, ncol = 5)
j = 1
for (i in 1:length(aresults)) {
  for (r in aresults[[i]]) {
    if (!is.null(r$dgp)) {
      mmses[j, 1] =  r$hs_mse / r$oracle_mse
      mmses[j+1, 1] = r$ss_mse / r$oracle_mse
      mmses[j+2, 1] = r$gl / r$oracle_mse
      mmses[j, 2] =  match(scenarios$p[r$dgp], levels(scenarios$p))
      mmses[j+1, 2] =  4
      mmses[j+2, 2] =  5
      mmses[j:(j+2), 3] =  match(scenarios$sparse[r$dgp], levels(scenarios$sparse))
      mmses[j:(j+2), 4] = scenarios$snr[r$dgp]
      mmses[j:(j+2), 5] = scenarios$corr[r$dgp]
      j = j + 3
    }
  }
}

df = data.frame(mmses)
colnames(df) = c("mse", "p", "sparse", "snr", "corr")
df = df[df$snr != 0, ]
df$snr = factor(df$snr , labels = c("'SNR = 5'", "'SNR = 20'"))
df$sparse = factor(df$sparse , labels = c("'high sparsity'", "'low sparsity'"))
df$p = factor(df$p , labels = c("strict", "optimal", "loose", "spikeSlabGAM","GAMSEL"))
df$corr = factor(df$corr , labels = c("rho*' = 0.0'", "rho*' = 0.7'"))

# without GAMSEL
df %>% filter(df$p != "GAMSEL") %>% 
  ggplot(aes(x = p, y = mse, group = p)) +
  geom_boxplot()  +
  scale_y_continuous(breaks = c(1,2)) + 
  coord_trans(y = "log10") + 
  xlab("method") + ylab(TeX('$\\log_{10}$(MSE / OracleMSE)'))+#"log(MSE / OracleMSE)") +
  facet_grid(sparse ~ snr + corr, labeller = label_parsed, scales = "free") +
  #ylim(0,2.5) + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

ggsave("../plots/sf-msea.pdf")
# with GAMSEL
df %>% ggplot(aes(x = p, y = mse, group = p)) + geom_boxplot()  +
  xlab("method") + ylab(TeX('$\\log_{10}$(MSE / OracleMSE)')) +
  scale_y_continuous(breaks = c(1,10,100,1000)) + 
  coord_trans(y = "log10") + 
  facet_grid(sparse ~ snr + corr, labeller = label_parsed, scales = "free") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))
ggsave("../plots/sf-mseb.pdf")

# MSE of function outputs

nsz = nrow(scenarios)
qs =  matrix(0, nrow = 2 * nsz * 100, ncol = 6)
j = 1
for (i in 1:length(aresults)) {
  for (r in aresults[[i]]) {
    if (!is.null(r$dgp)) {
      qs[j, 1] =  r$qnz
      qs[j + 1, 1] =  r$qz
      qs[j + 1, 6] =  0
      qs[j, 6] =  1
      qs[j:(j + 1), 2] =  match(scenarios$p[r$dgp], levels(scenarios$p))
      qs[j:(j + 1), 3] =  match(scenarios$sparse[r$dgp], levels(scenarios$sparse))
      qs[j:(j + 1), 4] = scenarios$snr[r$dgp]
      qs[j:(j+1), 5] = scenarios$corr[r$dgp]
      j = j + 2
    }
  }
}

df = drop_na(data.frame(qs))
colnames(df) = c("mse", "p", "sparse", "snr","corr", "zero")
df = df[df$p != 0, ]
df$snr = factor(df$snr , labels = c("'SNR = 5'", "'SNR = 20'"))
df$sparse = factor(df$sparse , labels = c("'high sparsity'", "'low sparsity'"))
df$p = factor(df$p , labels = c("strict", "optimal", "loose"))
df$corr = factor(df$corr , labels = c("rho*' = 0.0'", "rho*' = 0.7'"))
df$zero = factor(df$zero , labels = c("zero", "nonzero"))

ggplot(df, aes(x = p, y = mse)) +
  geom_boxplot(data = df, aes(fill = zero, color = zero)) +
  scale_fill_manual(name = "Output of function",
                    values = c("grey40", "lightskyblue")) +
  scale_color_manual(name = "Output of function", values = c("black", "blue")) +
  xlab("prior") + ylab(TeX("MSE of function output estimate $\\bar{|f - \\hat{f}|^2_2}$")) +
  facet_grid(sparse ~ snr + corr, labeller = label_parsed, scales = "free") +
  theme(legend.position = "bottom",
        legend.title=element_text(size=17),
        legend.text=element_text(size=16),
        axis.text.x=element_text(angle=45, hjust=1),
        axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))

ggsave("../plots/sf-quality.pdf", width=8)


