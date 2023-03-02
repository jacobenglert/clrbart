# Program Name: clrbart.R
# Author: Jacob Englert
# Date: 11 January 2023
# Purpose: Develop framework for performing conditional logistic regression
# using Bayesian Additive Regression Trees (BART) where only binary splits
# exist.

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(survival)
library(ars)
source('helper_fcns.R')

# Load Test Data ----------------------------------------------------------
data <- read_csv('Data/test.csv') %>%
  select(Y, Z, ID, starts_with('X'))

y <- data$Y
z <- data$Z
x <- as.matrix(data[,which(!(colnames(data) %in% c('ID','Y','Z')))])
strata <- data$ID

# Define clrbart function -------------------------------------------------
clrbart <- function(x, y, z, strata, sigma2.mu = 1, iter = 100,
                     init.tree = PLANT(x, y, z, sc1 = 0, strata, sigma2.mu), 
                     moves = c('BIRTH','DEATH','CHANGE'), 
                     move_probs = rep(1/3, 3),
                     min_leaf_size = 20){
  
  # Initialize tree
  tree <- init.tree$tree
  
  # Storage for proposed moves, accept/reject status, run time, and log-likelihood
  move_store <- character(iter)
  acc_store <- character(iter)
  time_store <- numeric(iter)
  logLik <- numeric(iter)
  mu_store <- matrix(nrow = iter, ncol = length(unique(strata)))
  
  # Begin MCMC
  mcmc.start <- Sys.time()
  for(i in 1:iter){
    start <- Sys.time()
    
    # Identify nodes in the current tree
    nodes <- which(lapply(tree, length) > 0)
    
    # If tree only consists of the root, force proposal to BIRTH
    if(length(nodes) == 1) move_probs <- c(1, 0, 0)
    
    # Randomly select a movement type
    move <- sample(moves, 1, prob = move_probs)
    
    # Depending on the movement type, generate a proposal and MH acceptance probability
    if(move == 'BIRTH'){
      prop_tree <- BIRTH(tree = tree, x = x, y = y, z = z, strata = strata, sigma2.mu = sigma2.mu, min_leaf_size = min_leaf_size)
      
      # If no BIRTH moves are possible then try another move
      if(min(prop_tree$new.nodes) == 0){
        move <- sample(moves, 1, prob = c(0, 0.5, 0.5))
        if(move == 'DEATH'){
          prop_tree <- DEATH(tree = tree, y = y, z = z, strata = strata, sigma2.mu = sigma2.mu)
          r <- rDEATH(tree, prop_tree$tree, pDEATH = 1/2, pBIRTH = 1/3)
        }
        if(move == 'CHANGE'){
          prop_tree <- CHANGE(tree = tree, x = x, y = y, z = z, strata = strata, sigma2.mu = sigma2.mu, min_leaf_size = min_leaf_size)
    
          # If no BIRTH or CHANGE moves are possible then force a death
          if(min(prop_tree$new.nodes) == 0){
            prop_tree <- DEATH(tree = tree, y = y, z = z, strata = strata, sigma2.mu = sigma2.mu)
            r <- rDEATH(tree, prop_tree$tree, pDEATH = 1, pBIRTH = 1/3)
          }
          r <- rCHANGE(tree, prop_tree$tree)
        }
      }
      r <- rBIRTH(tree, prop_tree$tree, pBIRTH = move_probs[1])
    }
    else if(move == 'DEATH'){
      prop_tree <- DEATH(tree = tree, y = y, z = z, strata = strata, sigma2.mu = sigma2.mu)
      r <- rDEATH(tree, prop_tree$tree)
    }
    else if(move == 'CHANGE'){
      prop_tree <- CHANGE(tree = tree, x = x, y = y, z = z, strata = strata, sigma2.mu = sigma2.mu, min_leaf_size = min_leaf_size)
      
      # If no CHANGE moves are possible then try another move
      if(min(prop_tree$new.nodes) == 0){
        move <- sample(moves, 1, prob = c(1/2, 1/2, 0))
        if(move == 'DEATH'){
          prop_tree <- DEATH(tree = tree, y = y, z = z, strata = strata, sigma2.mu = sigma2.mu)
          r <- rDEATH(tree, prop_tree$tree, pBIRTH = 1/2, pDEATH = 1/3)
        }
        if(move == 'BIRTH'){
          prop_tree <- BIRTH(tree = tree, x = x, y = y, z = z, strata = strata, sigma2.mu = sigma2.mu, min_leaf_size = min_leaf_size)
          
          # If no BIRTH or CHANGE moves are possible then force a death
          if(min(prop_tree$new.nodes) == 0){
            prop_tree <- DEATH(tree = tree, y = y, z = z, strata = strata, sigma2.mu = sigma2.mu)
            r <- rDEATH(tree, prop_tree$tree, pBIRTH = 1/3, pDEATH = 1)
          }
          r <- rBIRTH(tree, prop_tree$tree, pBIRTH = 1/2, pDEATH = 1/3)
        }
      }
      r <- rCHANGE(tree, prop_tree$tree)
    }
    
    # Accept or reject proposed tree via MH step
    if(log(runif(1)) < r){
      tree <- sampleM(tree = prop_tree$tree, nodes = prop_tree$new.nodes, y = y, z = z, strata = strata, sigma2.mu = sigma2.mu)
      sigma2.mu <- sampleS(tree = prop_tree$tree, sigma2.mu = sigma2.mu)
      acc <- 'Accepted'
      cat(move, 'move proposed and accepted.', prop_tree$message, '\n')
    } else{
      acc <- 'Rejected'
      cat(move, 'move proposed and rejected. \n')
    } 
    
    # Store movement, run time, acceptance, and log-likelihood information
    move_store[i] <- move
    time_store[i] <- Sys.time() - start
    acc_store[i] <- acc
    logLik[i] <- sum(unlist(lapply(tree, function(l) l$l*l$logLik)))
    
    mu_store[i,] <- getBARTfits(tree, strata)
    
    # Reset move probabilities
    move_probs <- c(1/3, 1/3, 1/3)
    
    # Update user on progress
    if((i %% 10) == 0) cat(paste0('Iteration ', i, ' completed. Total time elapsed: ', round(difftime(Sys.time(), mcmc.start, units = 'auto'), 1), ' seconds. \n'))
  }
  
  return(list(tree = tree, mu = mu_store, tree_detail = data.frame(move = move_store, time = time_store, acc = acc_store, logLik = logLik)))
}

# Sandbox -----------------------------------------------------------------

# Run tree
set.seed(42)
tree <- clrbart(x = x, y = y, z = z, strata = strata, iter = 200, sigma2.mu = .01)

# Summarize tree
cbind(id = which(lapply(tree$tree, length) > 0), do.call(rbind, tree$tree))
sum(unlist(lapply(tree$tree, function(l) l$l*l$logLik)))

# Display tree
data.frame(cbind(strata, x)) %>%
  distinct() %>%
  cbind(t(tree$mu)) %>%
  select(-strata) %>%
  distinct() %>%
  pivot_longer(cols = -starts_with('X'), names_to = 'iter', values_to = 'mu') %>%
  mutate(iter = as.numeric(iter)) %>%
  ggplot(aes(x = iter, y = mu, color = interaction(X1, X2))) +
  geom_line() +
  #geom_step(direction = 'hv') +
  theme_bw() +
  labs(x = 'MCMC Iteration',
       y = 'Log-OR')


# Visualize iteration time, likelihood, and moves
data.frame(tree$tree_detail) %>%
  mutate(iter = row_number()) %>%
  pivot_longer(cols = c(logLik, time), names_to = 'Metric', values_to = 'Value') %>%
  ggplot(aes(x = iter, y = Value)) +
  geom_line() +
  geom_point(aes(color = acc, shape = move), size = 2) +
  theme_bw() +
  facet_wrap(~Metric, ncol = 1, scales = 'free')


# Visualize the fisher scoring results
nodes <- which(!(tree$tree %in% list(NULL)))
mhat <- unlist(lapply(tree$tree, '[[', 'm'))
vhat <- unlist(lapply(tree$tree, '[[', 'v'))
mseq <- seq(min(mhat - 4*vhat) - 5, max(mhat + 4*vhat) + 5, 0.1)
lik <- matrix(ncol = length(nodes), nrow = length(mseq))
for(i in 1:length(nodes)){
  node <- nodes[i]
  id <- tree$tree[[node]]$DataIDs
  lik[,i] <- sapply(mseq, function(m) compLogLik(y[id], sc = z[id] * m, strata[id], max.win = tree$tree[[node]]$mw, na.locs = tree$tree[[node]]$nas))
}

lik %>%
  as.data.frame() %>%
  mutate(m = mseq) %>%
  pivot_longer(cols = -m, names_to = 'node', values_to = 'lik') %>%
  # group_by(node) %>%
  # mutate(lik = (lik - min(lik)) / (max(lik) - min(lik))) %>%
  # ungroup() %>%
  mutate(mhat = rep(mhat, times = length(mseq)),
         node = rep(nodes, times = length(mseq))) %>%
  ggplot(aes(x = m, y = lik, color = factor(node))) +
  geom_line(alpha = 0.4) +
  geom_vline(aes(xintercept = mhat, color = factor(node)), lty = 2) +
  theme_bw() +
  facet_wrap(~node, scales = 'free')
