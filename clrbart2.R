# Program Name: clrbart2.R
# Author: Jacob Englert
# Date: 28 February 2023
# Purpose: Develop framework for performing conditional logistic regression
# using Bayesian Additive Regression Trees (BART) where only binary splits
# exist. Extends former clrbart() functionality to include another set of
# fixed effects (binary or continuous) in the systematic component of the model.

# Load Packages -----------------------------------------------------------
packages <- c('tidyverse','survival','ars','progress')
install.packages(setdiff(packages, rownames(installed.packages())))
rm(packages)
library(tidyverse)
library(survival)
library(ars)
library(progress)

# Define clrbart2 function ------------------------------------------------
clrbart2 <- function(w, x, y, z, strata, beta.corr = solve(t(x) %*% x),
                     iter = 1e4, sigma2.beta = 0.01, acc.prob = 0.35, thin = 1, burnin = iter/2,
                     sigma2.mu = 1, sigma2.beta.update.freq = 100,
                     moves = c('BIRTH','DEATH','CHANGE'), 
                     move_probs = rep(1/3, 3),
                     min_leaf_size = 20){
  
  cat('Initializing Model Parameters... \n')
  
  # Parameters needed to accelerate computation of log-likelihood on the entire dataset
  windows <- table(strata)
  max.win <- max(windows)
  na.locs <- max.win * unique(strata)[which(windows != max.win)]
  
  # Number of strata and confounders
  n <- length(unique(strata))
  p <- ncol(x)
  
  # Set up storage for MCMC results
  beta <- matrix(0, ncol = p, nrow = iter)
  colnames(beta) <- paste0('beta', 1:p)
  
  
  mu.lookup <- list()
  # mu <- matrix(0, ncol = n, nrow = iter)
  # colnames(mu) <- paste0('mu', 1:n)
  
  
  sigma2.beta.store <- matrix(0, nrow = iter)
  sigma2.mu.store   <- matrix(0, nrow = iter)
  move.store        <- matrix(character(iter), nrow = iter)
  acc.store         <- matrix(0, nrow = iter)
  time.store        <- matrix(0, nrow = iter)
  logLik.store      <- matrix(0, nrow = iter)
  
  # Initialize beta and mu
  m0 <- clogit(y ~ x + z + strata(strata))
  beta[1,] <- coef(m0)[1:p]
  
  mu.lookup[[1]] <- matrix(coef(m0)[p+1])
  colnames(mu.lookup[[1]]) <- 'mu'
  # mu[1,] <- coef(m0)[p+1]
  
  sigma2.beta.store[1,] <- sigma2.beta
  sigma2.mu.store[1,] <- sigma2.mu
  
  # Initialize tree structure
  tree <- PLANT(x = w, y = y, z = z, sc1 = x %*% beta[1,], strata = strata, sigma2.mu = sigma2.mu)$tree
  
  logLik.store[1,] <- sum(unlist(lapply(tree, function(l) l$l * l$logLik)))
  
  # Begin tracking acceptance probability to automatically update sigma2.beta
  n.acc <- 0  # running acceptance counter
  n.acc2 <- 0 # current acceptance counter (updates with sigma2.beta)
  
  # Begin MCMC
  cat('Running MCMC... \n')
  pb <- progress_bar$new(
    format = "Running MCMC [:bar] :elapsedfull",
    total = 100, clear = FALSE, width = 60)
  pb$tick()
  mcmc.start <- Sys.time()
  for(k in 1:(iter-1)){
    start <- Sys.time()
    
    # Step 1: Update beta
    # Betas (proposal and current)
    beta.prop <- as.numeric(MASS::mvrnorm(1, beta[k,], sigma2.beta * beta.corr))
    beta.curr <- as.numeric(beta[k,])
    
    # Log-likelihoods (proposal and current)
    # sc2 <- rep.int(mu[1,], windows)
    ljby <- setdiff(names(which(apply(mu.lookup[[k]], 2, var) != 0)),'mu')
    sc2 <- dplyr::left_join(as.data.frame(w), as.data.frame(mu.lookup[[k]]), by = ljby)$mu
    ll.prop <- compLogLik(y = y, sc = x %*% beta.prop + sc2, strata = strata, max.win = max.win, na.locs = na.locs)
    ll.curr <- compLogLik(y = y, sc = x %*% beta.curr + sc2, strata = strata, max.win = max.win, na.locs = na.locs)
    
    # Log-priors (proposal and current)
    lprior.prop <- mvtnorm::dmvnorm(beta.prop, sigma = 1e10*diag(p), log = TRUE)
    lprior.curr <- mvtnorm::dmvnorm(beta.curr, sigma = 1e10*diag(p), log = TRUE)
    
    # Transition kernels (these "cancel out" when proposal distribution is symmetric)
    # lq.prop <- mvtnorm::dmvnorm(beta.prop, beta.curr, beta.Sigma, log = TRUE)
    # lq.curr <- mvtnorm::dmvnorm(beta.curr, beta.prop, beta.Sigma, log = TRUE)
    
    # Acceptance ratio
    a <- min(0, ll.prop + lprior.prop - ll.curr - lprior.curr) # + lq.curr - lq.prop
    
    if(log(runif(1)) < a) {
      beta[k+1,] <- beta.prop
      if((k+1) > burnin) n.acc <- n.acc + 1
      n.acc2 <- n.acc2 + 1
    }
    else beta[k+1,] <- beta.curr
    
    # Update sigma2.beta
    if((k+1) %% sigma2.beta.update.freq == 0){
      if(n.acc2/sigma2.beta.update.freq < 0.001) sigma2.beta <- sigma2.beta*0.1
      else if(n.acc2/sigma2.beta.update.freq < 0.05) sigma2.beta <- sigma2.beta*0.5
      else if(n.acc2/sigma2.beta.update.freq < 0.20) sigma2.beta <- sigma2.beta*0.9
      else if(n.acc2/sigma2.beta.update.freq > 0.50) sigma2.beta <- sigma2.beta*1.1
      else if(n.acc2/sigma2.beta.update.freq > 0.75) sigma2.beta <- sigma2.beta*2.0
      else if(n.acc2/sigma2.beta.update.freq > 0.95) sigma2.beta <- sigma2.beta*10.0
      n.acc2 <- 0
    }
    sigma2.beta.store[k+1,] <- sigma2.beta
    
    xb <- x %*% beta[k+1,]
    
    # Step 2: Update mu and the tree structure
    # Identify nodes in the current tree
    nodes <- which(lapply(tree, length) > 0)
    
    # If tree only consists of the root, force proposal to BIRTH
    if(length(nodes) == 1) move_probs <- c(1, 0, 0)
    
    # Randomly select a movement type
    move <- sample(moves, 1, prob = move_probs)
    
    # Depending on the movement type, generate a proposal and MH acceptance probability
    if(move == 'BIRTH'){
      prop_tree <- BIRTH(tree = tree, x = w, y = y, z = z, sc1 = xb, strata = strata, sigma2.mu = sigma2.mu, min_leaf_size = min_leaf_size)
      
      # If no BIRTH moves are possible then try another move
      if(min(prop_tree$new.nodes) == 0){
        move <- sample(moves, 1, prob = c(0, 0.5, 0.5))
        if(move == 'DEATH'){
          prop_tree <- DEATH(tree = tree, y = y, z = z, sc1 = xb, strata = strata, sigma2.mu = sigma2.mu)
          r <- rDEATH(tree, prop_tree$tree, pDEATH = 1/2, pBIRTH = 1/3)
        }
        if(move == 'CHANGE'){
          prop_tree <- CHANGE(tree = tree, x = w, y = y, z = z, sc1 = xb, strata = strata, sigma2.mu = sigma2.mu, min_leaf_size = min_leaf_size)
          
          # If no BIRTH or CHANGE moves are possible then force a death
          if(min(prop_tree$new.nodes) == 0){
            prop_tree <- DEATH(tree = tree, y = y, z = z, sc1 = xb, strata = strata, sigma2.mu = sigma2.mu)
            r <- rDEATH(tree, prop_tree$tree, pDEATH = 1, pBIRTH = 1/3)
          }
          r <- rCHANGE(tree, prop_tree$tree)
        }
      }
      r <- rBIRTH(tree, prop_tree$tree, pBIRTH = move_probs[1])
    }
    else if(move == 'DEATH'){
      prop_tree <- DEATH(tree = tree, y = y, z = z, sc1 = xb, strata = strata, sigma2.mu = sigma2.mu)
      r <- rDEATH(tree, prop_tree$tree)
    }
    else if(move == 'CHANGE'){
      prop_tree <- CHANGE(tree = tree, x = w, y = y, z = z, sc1 = xb, strata = strata, sigma2.mu = sigma2.mu, min_leaf_size = min_leaf_size)
      
      # If no CHANGE moves are possible then try another move
      if(min(prop_tree$new.nodes) == 0){
        move <- sample(moves, 1, prob = c(1/2, 1/2, 0))
        if(move == 'DEATH'){
          prop_tree <- DEATH(tree = tree, y = y, z = z, sc1 = xb, strata = strata, sigma2.mu = sigma2.mu)
          r <- rDEATH(tree, prop_tree$tree, pBIRTH = 1/2, pDEATH = 1/3)
        }
        if(move == 'BIRTH'){
          prop_tree <- BIRTH(tree = tree, x = w, y = y, z = z, sc1 = xb, strata = strata, sigma2.mu = sigma2.mu, min_leaf_size = min_leaf_size)
          
          # If no BIRTH or CHANGE moves are possible then force a death
          if(min(prop_tree$new.nodes) == 0){
            prop_tree <- DEATH(tree = tree, y = y, z = z, sc1 = xb, strata = strata, sigma2.mu = sigma2.mu)
            r <- rDEATH(tree, prop_tree$tree, pBIRTH = 1/3, pDEATH = 1)
          }
          r <- rBIRTH(tree, prop_tree$tree, pBIRTH = 1/2, pDEATH = 1/3)
        }
      }
      r <- rCHANGE(tree, prop_tree$tree)
    }
    
    # Accept or reject proposed tree via MH step
    if(log(runif(1)) < r){
      tree <- sampleM(tree = prop_tree$tree, nodes = prop_tree$new.nodes, y = y, z = z, sc1 = xb, strata = strata, sigma2.mu = sigma2.mu)
      sigma2.mu <- sampleS(tree = prop_tree$tree, sigma2.mu = sigma2.mu)
      acc <- 1
    } else acc <- 0
    
    # Update mu
    # mu[k+1,] <- getBARTfits(tree, strata)
    mu.lookup[[k+1]] <- getLookupTable(tree, w, strata)
    sigma2.mu.store[k+1,] <- sigma2.mu
    
    # Store movement, run time, acceptance, and log-likelihood information
    move.store[k+1,] <- move
    time.store[k+1,] <- Sys.time() - start
    acc.store[k+1,] <- acc

    # Identify which nodes currently exist in the tree
    nodes <- which(lapply(tree, length) > 0)
    
    # Identify which of the nodes are leaves
    leaves <- nodes[unlist(lapply(tree[nodes], '[[', 'l')) == 1]
    
    # Update log-likelihood
    logLik.store[k+1,] <- sum(unlist(lapply(tree[leaves], function(l){
      id <- l$DataIDs
      compLogLik(y = y[id], sc = xb[id] + l$mu * z[id], strata = strata[id], max.win = l$mw, na.locs = l$nas)
    })))
    
    # Reset move probabilities
    move_probs <- c(1/3, 1/3, 1/3)
    
    # Update user on progress
    # if(((k+1) %% 10) == 0) cat(paste0('Iteration ', k+1, ' completed. Total time elapsed: ', round(difftime(Sys.time(), mcmc.start, units = 'auto'), 1), ' seconds. \n'))
    pb$tick()
  }
  
  # Apply burn-in and thinning before returning result (this can be more efficient with list)
  cat('Preparing results... \n')
  post <- list(beta = beta, 
               #mu = mu, 
               sigma2.beta = sigma2.beta.store, sigma2.mu = sigma2.mu.store,
               move = move.store, time = time.store, acc = acc.store, logLik = logLik.store)
  if(burnin > 0) post <- lapply(post, function(x) x[(burnin+1):iter, ])
  if(thin > 1) post <- lapply(post, function(x) x[seq(1, iter, thin), ])
  
  if(burnin > 0) mu.lookup <- mu.lookup[(burnin+1):iter]
  if(thin > 1) mu.lookup <- mu.lookup[seq(1, iter, thin)]
  
  cat('Finished! \n')
  
  # return(list(post = post, mu.lookup = mu.lookup, tree = tree, acc.prob = n.acc/(iter - burnin)))
  return(list(post = post, tree = tree, acc.prob = n.acc/(iter - burnin)))
}
