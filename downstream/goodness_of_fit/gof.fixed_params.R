

library(tidyverse)
library(MASS)

SEED <- 4321
set.seed(SEED)

path_m_tot <- fs::path("/Users/ieo6983/Desktop/phylo/gof/DP.csv")
path_m_suc <- fs::path("/Users/ieo6983/Desktop/phylo/gof/AD.csv")


##


# Matrices cells x variants 
m_tot <- read.csv(path_m_tot, header = T, row.names = 1) %>% as.matrix(.)
m_suc <- read.csv(path_m_suc, header = T, row.names = 1) %>% as.matrix(.)
AF <- round(m_suc / m_tot,3)

#
AF_mean <- round(apply(AF, MARGIN = 2, mean),3)
suc_mean <- round(apply(m_suc, MARGIN = 2, mean), 3)
suc_var <- round(apply(m_suc, MARGIN = 2, var), 3)
suc_sd <- round(apply(m_suc, MARGIN = 2, sd), 3)
trials_mean <- round(apply(m_tot, MARGIN = 2, mean))


##


# Parametrize with avg(AF) or other
# - bin, nb, poi
# One dist x variants

par(mfrow = c(1,3))
for(i in 1:dim(m_suc)[2]){
  successes <- m_suc[,i]
  n_trials <- trials_mean[i]
  
  # QQ plots
  
  # Binomial - prob is avg. AF
  qqplot(qbinom(ppoints(length(successes)), size = n_trials, prob = AF_mean[i]), successes, 
         main = paste0("Binomial - var: ", i), ylab = "observed", xlab = "theoretical")
  abline(0, 1, col = "red")

  
  # Poisson - lambda is avg. number of successes
  qqplot(qpois(ppoints(length(successes)), lambda = suc_mean[i]), successes, 
         main = paste0("Poisson - var: ", i), ylab = "observed", xlab = "theoretical")
  abline(0, 1, col = "red")
  
  # Negative binomial - size is sd of successes; mu is avg. number of successes 
  # Any parametrization of nb using AF_mean[i] did not work for me
  qqplot(qnbinom(ppoints(length(successes)), size = suc_sd[i], mu = suc_mean[i]), successes, 
         main = paste0("Negative binomial - var: ", i), ylab = "observed", xlab = "theoretical")
  abline(0, 1, col = "red")


  # Chi-squared test
  
  #
  
  expected_counts_binom <- dbinom(0:(length(successes)-1), size = n_trials, prob = AF_mean[i]) 
  chisq_binom <- sum((successes - expected_counts_binom)^2 / expected_counts_binom, na.rm = T)
  p_value_binom <- pchisq(chisq_binom, df = length(successes) - 1, lower.tail = FALSE)  
  print(paste0("Chi-squared binomial: ", chisq_binom)) 
  print(paste0("P-value: ", p_value_binom))
  
  expected_counts_poi <- dpois(0:(length(successes)-1), lambda = suc_mean[i]) 
  chisq_poi <- sum((successes - expected_counts_poi)^2 / expected_counts_poi, na.rm = T)
  p_value_poi <- pchisq(chisq_poi, df = length(successes) - 1, lower.tail = FALSE)  
  print(paste0("Chi-squared poisson: ", chisq_poi)) 
  print(paste0("P-value: ", p_value_poi))
  
  expected_counts_nb <- dnbinom(0:(length(successes)-1), size = suc_sd[i], mu = suc_mean[i]) 
  chisq_nb <- sum((successes - expected_counts_nb)^2 / expected_counts_nb, na.rm = T)
  p_value_nb <- pchisq(chisq_nb, df = length(successes) - 1, lower.tail = FALSE)  
  print(paste0("Chi-squared negative binomial: ", chisq_nb)) 
  print(paste0("P-value: ", p_value_nb))
  
  #
  
}


##


# Parametrize with avg(AF) or other
# - bin, nb, poi
# One dist x cell

#
margin <- 1
AF_mean <- round(apply(AF, MARGIN = margin, mean),3)
suc_mean <- round(apply(m_suc, MARGIN = margin, mean), 3)
suc_var <- round(apply(m_suc, MARGIN = margin, var), 3)
suc_sd <- round(apply(m_suc, MARGIN = margin, sd), 3)
trials_mean <- round(apply(m_tot, MARGIN = margin, mean))

par(mfrow = c(1,3))
for(i in 1:dim(m_suc)[1]){
  successes <- m_suc[i,]
  n_trials <- trials_mean[i]
  
  # QQ plots
  
  # Binomial - prob is avg. AF
  qqplot(qbinom(ppoints(length(successes)), size = n_trials, prob = AF_mean[i]), successes, 
         main = paste0("Binomial - cell: ", i), ylab = "observed", xlab = "theoretical")
  abline(0, 1, col = "red")
  
  
  # Poisson - lambda is avg. number of successes
  qqplot(qpois(ppoints(length(successes)), lambda = suc_mean[i]), successes, 
         main = paste0("Poisson - cell: ", i), ylab = "observed", xlab = "theoretical")
  abline(0, 1, col = "red")
  
  # Negative binomial - size is sd of successes; mu is avg. number of successes 
  # Any parametrization of nb using AF_mean[i] did not work for me
  qqplot(qnbinom(ppoints(length(successes)), size = suc_sd[i], mu = suc_mean[i]), successes, 
         main = paste0("Negative binomial - cell: ", i), ylab = "observed", xlab = "theoretical")
  abline(0, 1, col = "red")
  
  
  # Chi-squared test
  
  #
  
  expected_counts_binom <- dbinom(0:(length(successes)-1), size = n_trials, prob = AF_mean[i]) 
  chisq_binom <- sum((successes - expected_counts_binom)^2 / expected_counts_binom, na.rm = T)
  p_value_binom <- pchisq(chisq_binom, df = length(successes) - 1, lower.tail = FALSE)  
  print(paste0("Chi-squared binomial: ", chisq_binom)) 
  print(paste0("P-value: ", p_value_binom))
  
  expected_counts_poi <- dpois(0:(length(successes)-1), lambda = suc_mean[i]) 
  chisq_poi <- sum((successes - expected_counts_poi)^2 / expected_counts_poi, na.rm = T)
  p_value_poi <- pchisq(chisq_poi, df = length(successes) - 1, lower.tail = FALSE)  
  print(paste0("Chi-squared poisson: ", chisq_poi)) 
  print(paste0("P-value: ", p_value_poi))
  
  expected_counts_nb <- dnbinom(0:(length(successes)-1), size = suc_sd[i], mu = suc_mean[i]) 
  chisq_nb <- sum((successes - expected_counts_nb)^2 / expected_counts_nb, na.rm = T)
  p_value_nb <- pchisq(chisq_nb, df = length(successes) - 1, lower.tail = FALSE)  
  print(paste0("Chi-squared negative binomial: ", chisq_nb)) 
  print(paste0("P-value: ", p_value_nb))
  
  #
  
}


##

