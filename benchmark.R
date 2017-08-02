##### Data Simulation #####
n <- 50 # number of trajectories per noise level
noises <- seq(0.2, 3, 0.2)

# Phase shifted
ps <- multi_sims(type = "ps", noises = noises, n = n)
# Phase shifted with trend
pst <- multi_sims(type = "pst", noises = noises, n = n, slope = 0.1)
# Amplitude noise
na <- multi_sims(type = "na", noises = noises, n = n)


##### Distance to Mean #####
cond <- "noise"
ti <- "Time"
mea <- "value"
lab <- "variable"

ps_mean <- dist_mean(data = ps, condition = cond, tcol = ti, measure = mea, label = lab, return.mean = F)
pst_mean <- dist_mean(data = pst, condition = cond, tcol = ti, measure = mea, label = lab, return.mean = F)
na_mean <- dist_mean(data = na, condition = cond, tcol = ti, measure = mea, label = lab, return.mean = F)


##### Pairwise Distances #####
# this will take about 8min
system.time({
  ps_pw <- all_pairwise_stats(data = ps, condition = cond, label = lab, measure = mea, k_roll_mean = 5)
  pst_pw <- all_pairwise_stats(data = pst, condition = cond, label = lab, measure = mea, k_roll_mean = 5)
  na_pw <- all_pairwise_stats(data = na, condition = cond, label = lab, measure = mea, k_roll_mean = 5)
})


##### Plots #####
library(ggplot2)
library(gridExtra)

theme_update(plot.title = element_text(hjust = 0.5), text = element_text(size=15))

# Euclidean to Mean
p1 <- ggplot(ps_mean, aes(x= noise, y = euclid_to_mean)) + geom_boxplot(aes(group = noise)) + ggtitle("Phase Shifted")
p2 <- ggplot(pst_mean, aes(x= noise, y = euclid_to_mean)) + geom_boxplot(aes(group = noise)) + ggtitle("Phase Shifted with Trend")
p3 <- ggplot(na_mean, aes(x= noise, y = euclid_to_mean)) + geom_boxplot(aes(group = noise)) + ggtitle("Noisy Amplitude")
grid.arrange(p1, p2, p3, ncol = 3)

# Overlap
p1 <- ggplot(ps_pw, aes(x = noise, y = Overlap)) + geom_boxplot(aes(group = noise)) + ggtitle("Phase Shifted") + scale_y_continuous(limits = c(0,1))
p2 <- ggplot(pst_pw, aes(x = noise, y = Overlap)) + geom_boxplot(aes(group = noise)) + ggtitle("Phase Shifted with Trend") + scale_y_continuous(limits = c(0,1))
p3 <- ggplot(na_pw, aes(x = noise, y = Overlap)) + geom_boxplot(aes(group = noise)) + ggtitle("Noisy Amplitude") + scale_y_continuous(limits = c(0,1))
grid.arrange(p1, p2, p3, ncol = 3)

# Pearson correlation
p1 <- ggplot(ps_pw, aes(x = noise, y = Pearson)) + geom_boxplot(aes(group = noise)) + ggtitle("Phase Shifted") + scale_y_continuous(limits = c(-1,1))
p2 <- ggplot(pst_pw, aes(x = noise, y = Pearson)) + geom_boxplot(aes(group = noise)) + ggtitle("Phase Shifted with Trend") + scale_y_continuous(limits = c(-1,1))
p3 <- ggplot(na_pw, aes(x = noise, y = Pearson)) + geom_boxplot(aes(group = noise)) + ggtitle("Noisy Amplitude") + scale_y_continuous(limits = c(-1,1))
grid.arrange(p1, p2, p3, ncol = 3)

# Spearman correlation
p1 <- ggplot(ps_pw, aes(x = noise, y = Spearman)) + geom_boxplot(aes(group = noise)) + ggtitle("Phase Shifted") + scale_y_continuous(limits = c(-1,1))
p2 <- ggplot(pst_pw, aes(x = noise, y = Spearman)) + geom_boxplot(aes(group = noise)) + ggtitle("Phase Shifted with Trend") + scale_y_continuous(limits = c(-1,1))
p3 <- ggplot(na_pw, aes(x = noise, y = Spearman)) + geom_boxplot(aes(group = noise)) + ggtitle("Noisy Amplitude") + scale_y_continuous(limits = c(-1,1))
grid.arrange(p1, p2, p3, ncol = 3)

# Kendall correlation
p1 <- ggplot(ps_pw, aes(x = noise, y = Kendall)) + geom_boxplot(aes(group = noise)) + ggtitle("Phase Shifted") + scale_y_continuous(limits = c(-1,1))
p2 <- ggplot(pst_pw, aes(x = noise, y = Kendall)) + geom_boxplot(aes(group = noise)) + ggtitle("Phase Shifted with Trend") + scale_y_continuous(limits = c(-1,1))
p3 <- ggplot(na_pw, aes(x = noise, y = Kendall)) + geom_boxplot(aes(group = noise)) + ggtitle("Noisy Amplitude") + scale_y_continuous(limits = c(-1,1))
grid.arrange(p1, p2, p3, ncol = 3)
