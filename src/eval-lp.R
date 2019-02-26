library("latex2exp")
library("tidyverse")
library("ggplot2")

source("util-lp.R")

load(file = get_lp_file())
scenarios = get_lp_scenarios()

aresults = list()
for (i in 1:nrow(scenarios)) {
  aresults[[i]] = list()
}

for (r in results) {
  aresults[[r$dgp]][[r$j]] = r
}

# MSE

nsz = nrow(scenarios)
mmses =  matrix(0, nrow = nsz * 100, ncol = 4)
j = 1
for (i in 1:length(aresults)) {
  for (r in aresults[[i]]) {
    if (!is.null(r$dgp)) {
      mmses[j, 1] =  r$ghs_mse / r$oracle_mse
      mmses[j, 2] =  match(scenarios$p[r$dgp], levels(scenarios$p))
      mmses[j, 3] =  match(scenarios$sparse[r$dgp], levels(scenarios$sparse))
      mmses[j, 4] = scenarios$snr[r$dgp]
      j = j + 1
    }
  }
}

df = drop_na(data.frame(mmses))
colnames(df) = c("mse", "p", "sparse", "snr")
df = df[df$p != 0,]
df$snr = factor(df$snr , labels = c("SNR = 0.1", "SNR = 1", "SNR = 5"))
df$sparse = factor(df$sparse , labels = c("high sparsity", "low sparsity"))
df$p = factor(df$p , labels = c("strict", "optimal", "loose"))

ggplot(df, aes(x = p, y = mse, group = p)) +
  geom_boxplot()  +
  xlab("prior") + ylab("MSE / OracleMSE") +
  facet_grid(sparse ~ snr, labeller = label_value) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12)
  )
ggsave("../plots/lp-mse.pdf")

# Parameter MSE

k = 1
nsz = nrow(scenarios)
qs =  matrix(0, nrow = 4 * 2 * nsz * 100, ncol = 6)
j = 1
for (i in 1:length(aresults)) {
  for (r in aresults[[i]]) {
    if (!is.null(r$dgp)) {
      for (k in 2 * (0:2)) {
        qs[j + k, 1] =  r$mse_levels[1, k / 2 + 1]
        qs[j + k + 1, 1] =  r$mse_levels[2, k / 2 + 1]
      }
      qs[j + 6, 1] =  r$mse_nonzero
      qs[j + 7, 1] =  r$mse_zero
      
      qs[j:(j + 7), 2] = match(scenarios$p[r$dgp], levels(scenarios$p))
      qs[j:(j + 7), 3] = match(scenarios$sparse[r$dgp], levels(scenarios$sparse))
      qs[j:(j + 7), 4] = scenarios$snr[r$dgp]
      qs[j:(j + 7), 5] = rep(c(0, 1), 4)
      qs[j:(j + 7), 6] = sort(rep(1:4, 2))
      j = j + 8
    }
  }
}

df = drop_na(data.frame(qs))
colnames(df) = c("mse", "p", "sparse", "snr", "zero", "gsize")
df = df[df$p != 0,]
snrs = c("SNR = 0.1", "SNR = 1", "SNR = 5")
df$snr = factor(df$snr , labels = snrs, levels = c("0.1", "1", "5"))
df$sparse = factor(df$sparse , labels = c("low sparsity", "high sparsity"))
df$p = factor(df$p , labels = c("strict", "optimal", "loose"))
df$zero = factor(df$zero , labels = c("zero", "non-zero"))
df$gsize = factor(df$gsize, labels = c("aggregated", "small", "medium", "large"))

for (k in 1:3) {
  df %>% filter(snr == snrs[k]) %>%
    ggplot(aes(x = p, y = mse)) +
    geom_boxplot(aes(fill = zero, color = zero)) +
    scale_fill_manual(name = "Influence of parameters",
                      values = c("grey40", "lightskyblue")) +
    scale_color_manual(name = "Influence of parameters", values = c("black", "blue")) +
    xlab("prior") + ylab(TeX(
      "MSE of parameter estimate $\\bar{|\\beta - \\hat{\\beta}|^2_2}$"
    )) +
    facet_grid(sparse ~ gsize, labeller = label_value) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 17),
      legend.text = element_text(size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 17),
      strip.text.x = element_text(size = 15),
      strip.text.y = element_text(size = 15)
    )
  ggsave(paste0("../plots/lp-paramMSE", as.character(k), ".pdf"),
         width = 8)
}
