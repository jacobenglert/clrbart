# Program Name: analysis.R
# Author: Jacob Englert
# Date: 1 March 2023
# Purpose: Run tests of the clrbart2 function defined in clrbart.R using data
# simulated from CCO_sim_confounders.R

# Load Packages -----------------------------------------------------------
library(tidyverse)

# Load clrbart2 function --------------------------------------------------
source('https://raw.githubusercontent.com/jacobenglert/clrbart/main/clrbart2.R')
source('https://raw.githubusercontent.com/jacobenglert/clrbart/main/clrbart_helpers.R')

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
myfit <- clrbart2(w, x, y, z, strata, iter = 25000, burnin = 5000)
saveRDS(myfit, file = "MCMCoutput/xwz_20000_10MAR2023.RData")


# Analyze fit -------------------------------------------------------------

myfit <- readRDS(file = "MCMCoutput/xwz_20000_10MAR2023.RData")

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

# Lookup table method
mu.df <- lapply(myfit$mu.lookup, as.data.frame) # Convert lookup matrices to data.frames
mu.post <- bind_rows(mu.df, .id = 'iter') # Bind the lookup data.frames together w/ iteration ID
mu.post[-1,] %>%
  mutate_at('iter', as.numeric) %>%
  ggplot(aes(x = iter, y = mu, color = interaction(W1, W2))) +
  geom_line() +
  theme_bw() +
  labs(x = 'MCMC Iteration', 
       y = 'Log-OR')

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


# Inference on the mu's ---------------------------------------------------

w0 <- data.frame(W1 = 1, W2 = 1, W3 = NA, W4 = NA, W5 = NA)
w1 <- data.frame(W1 = 1, W2 = 0, W3 = NA, W4 = NA, W5 = NA)
w_data <- as.data.frame(w)


w_compare <- function(w0, w1, w_data, mu, method = 'causal'){
  
  options(dplyr.summarise.inform = FALSE)
  
  # Names of effect modifiers for each group
  w0_vars <- names(w0)[!is.na(w0)]
  w1_vars <- names(w1)[!is.na(w1)]
  
  # All potential effect modifiers
  w_vars <- colnames(w_data)
  
  # Frequencies of combinations of effect modifiers in the data
  w_counts <- w_data %>%
    group_by(across(all_of(w_vars))) %>%
    summarise(n = n()) %>%
    ungroup()

  # If causal then g-comp, else subgroup averages
  if(method == 'causal'){
    w0_data <- w0[w0_vars] %>%
      cross_join(select(w_counts, -all_of(w0_vars))) %>%
      group_by(across(all_of(w_vars))) %>% 
      summarise(n = sum(n))
    
    w1_data <- w1[w1_vars] %>%
      cross_join(select(w_counts, -all_of(w1_vars))) %>%
      group_by(across(all_of(w_vars))) %>% 
      summarise(n = sum(n))
    
    w0_mu <- w0_data %>% 
      left_join(mu, by = w_vars, multiple = 'all') %>%
      group_by(iter) %>%
      summarise(mu = sum(n * mu) / sum(n)) %>%
      ungroup() %>%
      select(mu)
    
    w1_mu <- w1_data %>% 
      left_join(mu, by = w_vars, multiple = 'all') %>%
      group_by(iter) %>%
      summarise(mu = sum(n * mu) / sum(n)) %>%
      ungroup() %>%
      select(mu)
    
    CATE <- bind_rows(w0_mu = w0_mu, w1_mu = w1_mu, diff = w1_mu - w0_mu, .id = 'Group') 
    return(CATE)
  } else if(method == 'association'){
    
    w0_data <- w0[w0_vars] %>%
      left_join(w_counts, by = w0_vars, multiple = 'all') %>%
      group_by(across(all_of(w_vars))) %>% 
      summarise(n = sum(n))
    
    w1_data <- w1[w1_vars] %>%
      left_join(w_counts, by = w1_vars, multiple = 'all') %>%
      group_by(across(all_of(w_vars))) %>% 
      summarise(n = sum(n))
    
    w0_mu <- w0_data %>% 
      left_join(mu, by = w_vars, multiple = 'all') %>%
      group_by(iter) %>%
      summarise(mu = sum(n * mu) / sum(n)) %>%
      ungroup() %>%
      select(mu)
    
    w1_mu <- w1_data %>% 
      left_join(mu, by = w_vars, multiple = 'all') %>%
      group_by(iter) %>%
      summarise(mu = sum(n * mu) / sum(n)) %>%
      ungroup() %>%
      select(mu)
    
    EFF <- bind_rows(w0_mu = w0_mu, w1_mu = w1_mu, diff = w1_mu - w0_mu, .id = 'Group') 
    return(EFF)
  }
}

w0_w1_c <- w_compare(w0, w1, w_data, mu.post, method = 'causal')
w0_w1_a <- w_compare(w0, w1, w_data, mu.post, method = 'association')

summary(filter(w0_w1_c, Group == 'diff')$mu)
summary(filter(w0_w1_a, Group == 'diff')$mu)

# Plot
w0_w1_c %>%
  ggplot(aes(x = mu, fill = Group)) +
  geom_density(alpha = 0.7, color = 'black') +
  facet_wrap(~Group, ncol = 1)
ggsave(filename = 'Figures/CATEplot.png')

w0_w1_a %>%
  ggplot(aes(x = mu, fill = Group)) +
  geom_density(alpha = 0.7, color = 'black') +
  facet_wrap(~Group, ncol = 1)
ggsave(filename = 'Figures/Diffplot.png')


# Complete storage method
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