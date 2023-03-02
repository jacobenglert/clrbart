# Program Name: test_clrbart2.R
# Author: Jacob Englert
# Date: 1 March 2023
# Purpose: Run tests of the clrbart2 function defined in clrbart.R using data
# simulated from CCO_sim_confounders.R

# Load Packages -----------------------------------------------------------
library(dplyr)

# Load clrbart2 function --------------------------------------------------
source('clrbart2.R')

# Load Test Data ----------------------------------------------------------
data <- read_csv('Data/XWZ.csv', show_col_types = FALSE) %>%
  select(Y, Z, ID, starts_with('X'), starts_with('W'))

w <- as.matrix(select(data, starts_with('W')))
x <- as.matrix(select(data, starts_with('X')))
y <- data$Y
z <- data$Z
strata <- data$ID

# Fit Model ---------------------------------------------------------------
set.seed(2187)
myfit <- clrbart2(w, x, y, z, strata, iter = 10000, burnin = 0)
saveRDS(myfit, file = "MCMCMoutput/xwz_20000_01MAR2023.RData")

# Summarize tree
cbind(id = which(lapply(myfit$tree, length) > 0), do.call(rbind, myfit$tree))
sum(unlist(lapply(myfit$tree, function(l) l$l*l$logLik)))

# Visualize ---------------------------------------------------------------

# Betas
myfit$post$beta %>%
  as.data.frame() %>%
  mutate(iter = row_number()) %>%
  pivot_longer(cols = -iter, names_to = 'parm', values_to = 'value') %>%
  ggplot(aes(x = iter, y = value, color = parm)) +
  geom_line() +
  facet_wrap(~parm, scales = 'free', ncol = 1)

# Mus
mu.df <- lapply(myfit$mu.lookup, as.data.frame) # Convert lookup matrices to data.frames
mu.post <- bind_rows(mu.df, .id = 'iter') # Bind the lookup data.frames together w/ iteration ID
mu.post[-1,] %>%
  mutate_at('iter', as.numeric) %>%
  ggplot(aes(x = iter, y = mu, color = interaction(W1, W2))) +
  geom_line() +
  theme_bw() +
  labs(x = 'MCMC Iteration', 
       y = 'Log-OR')

# data.frame(cbind(strata, w)) %>%
#   distinct() %>%
#   cbind(t(myfit$post$mu)) %>%
#   select(-strata) %>%
#   distinct() %>%
#   pivot_longer(cols = -starts_with('W'), names_to = 'iter', values_to = 'mu') %>%
#   mutate(iter = as.numeric(iter)) %>%
#   ggplot(aes(x = iter, y = mu, color = interaction(W1, W2, W3))) +
#   geom_line() +
#   #geom_step(direction = 'hv') +
#   theme_bw() +
#   labs(x = 'MCMC Iteration',
#        y = 'Log-OR')

# Visualize iteration time, likelihood, and moves
data.frame(move = myfit$post$move,
           time = myfit$post$time,
           acc = myfit$post$acc,
           logLik = myfit$post$logLik) %>%
  mutate(iter = row_number()) %>%
  pivot_longer(cols = c(logLik, time), names_to = 'Metric', values_to = 'Value') %>%
  filter(iter %% 50 == 0) %>%
  ggplot(aes(x = iter, y = Value)) +
  geom_line() +
  geom_point(aes(color = factor(acc), shape = move), size = 2) +
  theme_bw() +
  facet_wrap(~Metric, ncol = 1, scales = 'free')

mean(myfit$post$acc)
