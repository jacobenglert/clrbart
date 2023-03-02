# Collection of helper functions for clrbart


# Define Helper Functions -------------------------------------------------

# getParents (obtain parent node id given child id)
getParent <- function(id){
  if(id == 1) return(numeric(0))
  if(id %% 2 == 0) p.id <- id / 2
  else p.id <- (id - 1) / 2
  return(p.id)
}

# getAncestors (obtain all ancestor ids given node id)
getAncestors <- function(id){
  a.id <- numeric()
  if(id == 1) return(numeric(0))
  while(length(id) > 0){
    id <- getParent(id)
    a.id <- c(a.id, id)
  }
  return(rev(a.id))
}

# getChildren (obtain left and right child id given parent id)
getChildren <- function(id){
  return(c(id*2, id*2 + 1))
}

# getSibling (obtain sibling id given node id)
getSibling <- function(id){
  c.ids <- getChildren(getParent(id))
  return(c.ids[which(c.ids != id)])
}

# resample (better version of default base R sample function)
resample <- function(x, ...) x[sample.int(length(x), ...)]

# getDataIDs (return a vector containing the row ids of observations for a given node)
getDataIDs <- function(node, tree, x){
  if(node == 1) return(1:nrow(x))
  else{
    ancestors <- getAncestors(node)
    a.rules <- paste(unlist(lapply(tree[c(ancestors, node)], function(x) x$rule)), collapse = ' & ')
    x %>%
      as.data.frame() %>%
      mutate(id = row_number()) %>% 
      filter(eval(str2lang(a.rules))) %>%
      select(id) %>%
      as_vector()
  }
}

# pSplit (return prior splitting probability)
pSplit <- function(d, gamma = 0.95, beta = 2) gamma * ((1 + d)^(-beta))

# getMWNA
getMWNA <- function(strata){
  strata <- as.numeric(as.factor(strata))
  windows <- table(strata)
  max.win <- max(windows)
  na.locs <- max.win * unique(strata)[which(windows != max.win)]
  return(list(mw = max.win, nas = na.locs))
}

# getBARTfits (returns fitted values for each unique strata)
getBARTfits <- function(tree, strata){
  fits <- numeric(length(unique(strata)))
  # Identify which nodes currently exist in the tree
  nodes <- which(lapply(tree, length) > 0)
  
  # Identify which of the nodes are leaves
  leaves <- nodes[unlist(lapply(tree[nodes], '[[', 'l')) == 1]
  
  for(l in leaves){
    id <- tree[[l]]$DataIDs
    fits[unique(strata[id])] <- tree[[l]]$mu
  }  
  return(fits)
}

# getBARTfits2 (returns fitted values for each unique strata)
getLookupTable <- function(tree, x, strata){

  # Identify which nodes currently exist in the tree
  nodes <- which(lapply(tree, length) > 0)
  
  # Create lookup table
  mu.lookup <- list()
  # if(length(nodes) == 1){
  #   mu.lookup[[1]] <- matrix(tree[[1]]$mu)
  #   colnames(mu.lookup[[1]]) <- 'mu'
  # } else{
    
    # Identify which of the nodes are leaves
    leaves <- nodes[unlist(lapply(tree[nodes], '[[', 'l')) == 1]
  
    # Identify variables that have been split on
    x_splits <- names(which(sapply(colnames(x), function(x) length(grep(x, unlist(lapply(tree, '[[', 'rule')))) > 0)))
  
    for(i in 1:length(leaves)){
      l <- leaves[i]
      id <- tree[[l]]$DataIDs
      xs <- as.matrix(dplyr::distinct(as.data.frame(cbind(x[id,], mu = tree[[l]]$mu)), across(all_of(colnames(x))), .keep_all = TRUE))
      colnames(xs) <- c(colnames(x), 'mu')
      mu.lookup[[i]] <- xs
    }  
  return(do.call(rbind, mu.lookup))
}

# Define Proposal Functions -----------------------------------------------

# Compute the CLR data log-likelihood for a given node
compLogLik <- function(y, sc, strata, max.win, na.locs){
  
  n.obs <- length(unique(strata))
  denoms <- NA * numeric(max.win*n.obs)
  
  denoms[-na.locs] <- sc
  tmp_d <- matrix(denoms, ncol = max.win, byrow = TRUE)
  
  logLik <- (denoms[-na.locs])[y == 1] - matrixStats::rowLogSumExps(tmp_d, na.rm = TRUE)
  return(sum(as.numeric(logLik)))
}

# Compute the gradient of the CLR likelihood for a given node
compGrd <- function(y, z, sc, strata, max.win, na.locs){
  
  n.obs <- length(unique(strata))
  nums <- denoms <- NA * numeric(max.win*n.obs)
  
  nums[-na.locs] <- z * exp(sc) # z * exp(z * mu)
  denoms[-na.locs] <- exp(sc) # exp(z * mu)
  
  tmp_n <- matrix(nums, ncol = max.win, byrow = TRUE)
  tmp_d <- matrix(denoms, ncol = max.win, byrow = TRUE)
  
  denom <- rowSums(tmp_d, na.rm = TRUE)
  num <- rowSums(tmp_n, na.rm = TRUE)
  
  grd <- z[y == 1] - num / denom
  
  return(sum(grd))
}

# Compute Fisher's information of the CLR likelihood for a given node
compFisher <- function(y, z, sc, strata, max.win, na.locs){
  
  n.obs <- length(unique(strata))
  a <- b <- c <- NA * numeric(max.win*n.obs)
  
  a[-na.locs] <- exp(sc)  # exp(z * mu)
  b[-na.locs] <- z * exp(sc)    # z * exp(z * mu)
  c[-na.locs] <- z^2 * exp(sc)  # z^2 * exp(z * mu)
  
  tmp_a <- matrix(a, ncol = max.win, byrow = TRUE)
  tmp_b <- matrix(b, ncol = max.win, byrow = TRUE)
  tmp_c <- matrix(c, ncol = max.win, byrow = TRUE)
  
  a <- rowSums(tmp_a, na.rm = TRUE)
  b <- rowSums(tmp_b, na.rm = TRUE)
  c <- rowSums(tmp_c, na.rm = TRUE)
  
  I <- ((a * c) - (b^2)) / (a^2)
  
  return(sum(I))
}

# mvCompute (compute the mean and variance for the proposal distribution)
mvCompute <- function(node, tree, y, z, sc1 = 0, strata, sigma2.mu, move, max.win, na.locs){
  
  # Initialize m
  if(move == 'BIRTH') m <- tree[[getParent(node)]]$m
  else if(move == 'DEATH'){
    # mean(unlist(lapply(tree[getChildren(node)], '[[', 'm')))
    ms <- unlist(lapply(tree[getChildren(node)], '[[', 'm'))
    ns <- unlist(lapply(tree[getChildren(node)], '[[', 'n'))
    m <- sum(ms * ns) / sum(ns)
  } 
  else if(move == 'PLANT') m <- coef(clogit(y ~ z + strata(strata)))
  
  # Use Fisher scoring to update m and v
  grd <- compGrd(y = y, z = z, sc = sc1 + z * m, strata = strata, max.win = max.win, na.locs = na.locs) - (m / sigma2.mu)
  I <- compFisher(y = y, z = z, sc = sc1 + z * m, strata = strata, max.win = max.win, na.locs = na.locs) + (1 / sigma2.mu)
  while(abs(grd) > sqrt(I)/10){
    m <- m + (grd / I)
    grd <- compGrd(y, z, sc1 + z * m, strata, max.win, na.locs) - (m / sigma2.mu)
    I <- compFisher(y, z, sc1 + z * m, strata, max.win, na.locs)  + (1 / sigma2.mu)
  }
  
  # Return the mean and standard deviation of the proposal distribution
  return(list(m = m, v = I^(-1/2)))
}

# sampleM (sample M from its full conditional using slice sampling)
sampleM <- function(tree, nodes, y, z, sc1 = numeric(length(y)), strata, sigma2.mu){
  
  # Log-likelihood and gradient functions needed for ars
  ll <- function(m, y, z, sc1, strata, sigma2.mu, max.win, na.locs){
    compLogLik(y, sc = sc1 + z * m, strata, max.win, na.locs) - ((m^2) / (2*sigma2.mu))
  }
  grd <- function(m, y, z, sc1, strata, sigma2.mu, max.win, na.locs){
    compGrd(y, z, sc1 + z * m, strata, max.win, na.locs) - (m / sigma2.mu)
  }
  
  ll.V <- Vectorize(FUN = ll, vectorize.args = 'm')
  grd.V <- Vectorize(FUN = grd, vectorize.args = 'm')
  
  for(n in nodes){
    id <- tree[[n]]$DataIDs
    tree[[n]]$mu <- ars::ars(1, ll.V, grd.V, y = y[id], z = z[id], sc1 = sc1[id], strata = strata[id], sigma2.mu = sigma2.mu, max.win = tree[[n]]$mw, na.locs = tree[[n]]$nas)
    tree[[n]]$logLik <- compLogLik(y = y[id], sc = sc1[id] + z[id] * tree[[n]]$mu, strata = strata[id], max.win = tree[[n]]$mw, na.locs = tree[[n]]$nas)
  }
  return(tree)
}

# SampleS (update sigma2.mu using its full conditional)
sampleS <- function(tree, sigma2.mu, alpha = 0.0001, beta = 0.0001){
  
  # Identify which nodes currently exist in the tree
  nodes <- which(lapply(tree, length) > 0)
  
  # Identify which of the nodes are leaves
  leaves <- nodes[unlist(lapply(tree[nodes], '[[', 'l')) == 1]
  
  n_mu <- length(leaves)
  mu <- unlist(lapply(tree[leaves], '[[', 'mu'))
  sigma2.mu <- 1 / rgamma(1, shape = n_mu / 2 + alpha, rate = sum((mu - 0)^2 / 2 + beta))
  
  return(sigma2.mu)
}


# Define BIRTH, DEATH, and CHANGE functions -------------------------------

# PLANT (initialize a new tree with a root node)
# BIRTH (create new leaf nodes for a current leaf node)
# DEATH (remove the leaf nodes from a non-grandparent branch node)
# CHANGE (change the splitting rule/variable for a non-grandparent branch node)
PLANT <- function(x, y, z, sc1 = numeric(length(y)), strata, sigma2.mu){
  
  # Instantiate tree
  tree <- list()
  
  # Instantiate root node (depth 0, not a branch, no rules, is a leaf, not a grandparent)
  tree[[1]] <- list(d = 0, b = 0, rule = character(0), l = 1, NOG = 0)
  
  # Every record belongs to the root node
  tree[[1]]$DataIDs <- id <- getDataIDs(1, tree, x)
  tree[[1]]$n <- length(unique(strata[id]))
  
  # Identify the maximum window size and NA locations for likelihood calculations
  mwna <- getMWNA(strata[id])
  tree[[1]]$mw <- mwna$mw
  tree[[1]]$nas <- mwna$nas
  
  # Compute m and v for the root node
  mv <- mvCompute(node = 1, tree = tree, y = y, z = z, sc1 = sc1, strata = strata, sigma2.mu = sigma2.mu, move = 'PLANT', max.win = mwna$mw, na.locs = mwna$nas)
  tree[[1]]$m <- mv$m
  tree[[1]]$v <- mv$v
  
  # Sample mu for the root node and update the log-likelihood
  tree[[1]]$mu <- rnorm(1, mv$m, mv$v)
  tree[[1]]$logLik <- compLogLik(y, sc = sc1 + z * tree[[1]]$mu, strata, max.win = mwna$mw, na.locs = mwna$nas)
  
  return(list(tree = tree, new.nodes = 1))
}
BIRTH <- function(tree, x, y, z, sc1 = numeric(length(y)), strata, sigma2.mu, min_leaf_size = 20, predictors = colnames(x)){
  
  # Identify which nodes currently exist in the tree
  nodes <- which(lapply(tree, length) > 0)
  
  # Identify which of the nodes are leaves
  leaves <- nodes[unlist(lapply(tree[nodes], '[[', 'l')) == 1]
  
  # Randomly select a leaf and check for valid splits. Randomly select from the
  # valid splits. If no valid split is found, repeat the process with another
  # leaf. If no leaves can be split, return the original tree with a note.
  no_split <- 1
  while(no_split == 1){
    
    if(length(leaves) == 0){
      message <- 'No valid splits were found. The tree has not been updated.'
      return(list(tree = tree, new.nodes = 0, message = message))
    }
    
    # Randomly select a leaf for splitting
    l <- resample(leaves, 1)
    
    # Identify which rules (variables) have already been used
    ancestors <- getAncestors(l)
    a.rules <- paste(unlist(lapply(tree[c(ancestors, l)], function(x) x$rule)), collapse = ' & ')
    used_x <- predictors[which(sapply(predictors, grep, a.rules) == 1)]
    
    # Identify new rules (variables) that would result in valid splits
    unused_x <- setdiff(predictors, used_x)
    
    if(l == 1) l_x <- as.data.frame(cbind(strata, x))
    else l_x <- dplyr::filter(as.data.frame(cbind(strata, x)), eval(str2lang(a.rules)))
    l_x_distinct <- as.matrix(dplyr::distinct(l_x)[, unused_x])
    new_node_sizes <- matrix(apply(l_x_distinct, 2, function(x) table(factor(x, levels = 0:1))), ncol = length(unused_x))
    colnames(new_node_sizes) <- unused_x
    
    # removes X's that would result in empty cell from consideration (this is a dimension check on the table)
    valid_x <- names(which(colSums(new_node_sizes > min_leaf_size) >= nrow(new_node_sizes))) 
    
    # Check another leaf if there are no available splits for this leaf
    if(length(valid_x) == 0){
      leaves <- setdiff(leaves, l)
      next
    }
    no_split <- 0
  }
  
  # Update rules for new split
  x_split <- sample(valid_x, 1)
  new_rule <- c(paste(x_split, '==', 0), paste(x_split, '==', 1))
  
  # Update original leaf node
  tree[[l]]$b <- 1
  tree[[l]]$l <- 0
  tree[[l]]$NOG <- 1
  
  # Update parent of original leaf node
  if(l != 1) tree[[getParent(l)]]$NOG <- 0
  
  # Create children of original leaf node
  l_children <- getChildren(l)
  for(i in 1:2){
    
    # Instantiate the child node
    c <- l_children[i]
    tree[[c]] <- list(d = tree[[l]]$d + 1, b = 0, rule = new_rule[i], l = 1, NOG = 0)
    
    # Identify records which belong to the child node
    tree[[c]]$DataIDs <- id <- getDataIDs(c, tree, x)
    tree[[c]]$n <- length(unique(strata[id]))
    
    # Identify the maximum window size and NA locations for likelihood calculations
    mwna <- getMWNA(strata[id])
    tree[[c]]$mw <- mwna$mw
    tree[[c]]$nas <- mwna$nas
    
    # Compute m and v for the child node
    mv <- mvCompute(node = c, tree = tree, y = y[id], z = z[id], sc1 = sc1[id], strata = strata[id], sigma2.mu = sigma2.mu, move = 'BIRTH', max.win = mwna$mw, na.locs = mwna$nas)
    tree[[c]]$m <- mv$m
    tree[[c]]$v <- mv$v
    
    # Sample mu for the child node and update the log-likelihood
    tree[[c]]$mu <- rnorm(1, mv$m, mv$v)
    tree[[c]]$logLik <- compLogLik(y = y[id], sc = sc1[id] + z[id] * tree[[c]]$mu, strata = strata[id], max.win = mwna$mw, na.locs = mwna$nas)
  }
  
  message <- paste0('In the proposed tree node ', l, ' has been split into node ', l_children[1], ' if ', new_rule[1], ' and node ', l_children[2], ' if ', new_rule[2], '.')
  return(list(tree = tree, new.nodes = l_children, message = message))
}
DEATH <- function(tree, y, z, sc1 = numeric(length(y)), strata, sigma2.mu){
  
  # Identify which nodes currently exist in the tree
  nodes <- which(lapply(tree, length) > 0)
  
  if(length(nodes) == 1){
    message <- 'Tree cannot be pruned since it only consists of the root node.'
    return(list(tree = tree, new.nodes = 0, message = message))
  }
  
  # Identify which of the nodes are NOG
  nogs <- nodes[unlist(lapply(tree[nodes], '[[', 'NOG')) == 1]
  
  # Randomly select a NOG node to prune back to
  b <- resample(nogs, 1)
  id <- tree[[b]]$DataIDs
  
  # Update original NOG node
  tree[[b]]$b <- 0
  tree[[b]]$l <- 1
  tree[[b]]$NOG <- 0
  
  # Recompute m and v for the NOG node (sigma2.mu may have changed since last time)
  mv <- mvCompute(node = b, tree = tree, y = y[id], z = z[id], sc1 = sc1[id], strata = strata[id], sigma2.mu = sigma2.mu, move = 'DEATH', max.win = tree[[b]]$mw, na.locs = tree[[b]]$nas)
  tree[[b]]$m <- mv$m
  tree[[b]]$v <- mv$v
  
  # Sample mu for the NOG node and update the log-likelihood
  tree[[b]]$mu <- rnorm(1, tree[[b]]$m, tree[[b]]$v)
  tree[[b]]$logLik <- compLogLik(y = y[id], sc = sc1[id] + z[id] * tree[[b]]$mu, strata = strata[id], max.win = tree[[b]]$mw, na.locs = tree[[b]]$nas)
  
  # Update parent of original NOG node (only if the DEATH makes the parent a NOG)
  if((b != 1) && (tree[[getSibling(b)]]$l == 1)) tree[[getParent(b)]]$NOG <- 1
  
  # Delete children of original NOG node
  b_children <- getChildren(b)
  tree[b_children] <- list(NULL)
  
  message <- paste0('In the proposed tree nodes ', b_children[1], ' and ', b_children[2], ' have been pruned back to node ', b, '.')
  
  return(list(tree = tree, new.nodes = b, message = message))
  
}
CHANGE <- function(tree, x, y, z, sc1 = numeric(length(y)), strata, sigma2.mu, min_leaf_size = 20, predictors = colnames(x)){
  
  # Identify which nodes currently exist in the tree
  nodes <- which(lapply(tree, length) > 0)
  
  if(length(nodes) == 1){
    message <- 'Tree cannot be changed since it only consists of the root node.'
    cat(message)
    return(list(tree = tree, new.nodes = 0, message = message))
  }
  
  # Identify which of the nodes are NOG
  nogs <- nodes[unlist(lapply(tree[nodes], '[[', 'NOG')) == 1]
  
  # Randomly select a NOG and check for valid changes. Randomly select from the
  # valid changes. If no valid change is found, repeat the process with another
  # NOG. If no NOGS can be changed, return the original tree with a note.
  no_change <- 1
  while(no_change == 1){
    
    if(length(nogs) == 0){
      message <- 'No valid changes were found. The tree has not been updated.'
      cat(message)
      return(list(tree = tree, new.nodes = 0, message = message))
    }
    
    # Randomly select a NOG node to change
    b <- resample(nogs, 1)
    
    # Identify which variables have already been used
    ancestors <- getAncestors(b)
    a.rules <- paste(unlist(lapply(tree[c(ancestors, b)], function(x) x$rule)), collapse = ' & ')
    used_x_above <- predictors[which(sapply(predictors, grep, a.rules) == 1)]
    
    lL <- getChildren(b)[1]
    used_x_below <- predictors[which(sapply(predictors, grep, tree[[lL]]$rule) == 1)]
    
    # Identify new rules (variables) that would result in valid splits
    unused_x <- setdiff(predictors, c(used_x_above, used_x_below))
    if(b == 1) b_x <- as.data.frame(cbind(strata, x))
    else b_x <- dplyr::filter(as.data.frame(cbind(strata, x)), eval(str2lang(a.rules)))
    b_x_distinct <- as.matrix(dplyr::distinct(b_x)[, unused_x])
    new_node_sizes <- matrix(apply(b_x_distinct, 2, function(x) table(factor(x, levels = 0:1))), nrow = 2, ncol = length(unused_x))
    colnames(new_node_sizes) <- unused_x
    
    # removes X's that would result in empty cell from consideration (this is a dimension check on the table)
    valid_x <- names(which(colSums(new_node_sizes > min_leaf_size) == nrow(new_node_sizes))) 
    
    # Check another NOG if there are no available changes for this NOG
    if(length(valid_x) == 0){
      nogs <- setdiff(nogs, b)
      next
    }
    no_change <- 0
  }
  
  # Update rules for new split
  x_change <- sample(valid_x, 1)
  new_rule <- c(paste(x_change, '==', 0), paste(x_change, '==', 1))
  
  # Update the leaf nodes off of the original NOG node b
  b_children <- getChildren(b)
  
  for(i in 1:2){
    
    # Update the rule for the child node
    c <- b_children[i]
    tree[[c]]$rule <- new_rule[i]
    
    # Identify records which belong to the child node
    tree[[c]]$DataIDs <- id <- getDataIDs(c, tree, x)
    tree[[c]]$n <- length(unique(strata[id]))
    
    # Identify the maximum window size and NA locations for likelihood calculations
    mwna <- getMWNA(strata[id])
    tree[[c]]$mw <- mwna$mw
    tree[[c]]$nas <- mwna$nas
    
    # Compute m and v for the child node
    mv <- mvCompute(node = c, tree = tree, y = y[id], z = z[id], sc1 = sc1[id], strata = strata[id], sigma2.mu = sigma2.mu, move = 'BIRTH', max.win = mwna$mw, na.locs = mwna$nas)
    tree[[c]]$m <- mv$m
    tree[[c]]$v <- mv$v
    
    # Sample mu for the child node and update the log-likelihood
    tree[[c]]$mu <- rnorm(1, mv$m, mv$v)
    tree[[c]]$logLik <- compLogLik(y = y[id], sc = sc1[id] + z[id] * tree[[c]]$mu, strata = strata[id], max.win = mwna$mw, na.locs = mwna$nas)
  }
  
  message <- paste0('In the proposed tree the splitting rule for node ', b, ' into nodes ', b_children[1], ' and ', b_children[2], ' is now based on ', x_change, '.')
  
  return(list(tree = tree, new.nodes = b_children, message = message))
}

# Define Metropolis-Hastings Acceptance Ratio Functions -------------------
rBIRTH <- function(curr_tree, prop_tree, pBIRTH = 1/3, pDEATH = 1/3){
  
  # Identify the proposed leaves and the parent node
  curr_leaves <- which(lapply(curr_tree, length) > 0)
  prop_leaves <- which(lapply(prop_tree, length) > 0)
  new_leaves <- setdiff(prop_leaves, curr_leaves)
  lL <- new_leaves[1]; lR <- new_leaves[2]; l <- getParent(lL)
  
  # Retrieve the log-likelihood for the proposed leaves and the parent node
  prop_logLik <- sum(vapply(c(lL, lR), function(l) prop_tree[[l]]$logLik, FUN.VALUE = 1))
  curr_logLik <- curr_tree[[l]]$logLik
  
  # Compute the splitting prior
  prop_rho <- pSplit(d = prop_tree[[lL]]$d)
  curr_rho <- pSplit(d = prop_tree[[l]]$d)
  
  # Compute the log-likelihood of the entire proposal
  gDEATH <- dnorm(curr_tree[[l]]$mu, curr_tree[[l]]$m, curr_tree[[l]]$v, log = TRUE)
  gBIRTH <- sum(unlist(lapply(prop_tree[c(lL, lR)], function(l) dnorm(l$mu, l$m, l$v, log = TRUE))))
  
  # Calculate the M-H acceptance ratio on the log scale
  a1 <- log(curr_rho) + 2 * log(1 - prop_rho) - log(1 - curr_rho)
  a2 <- prop_logLik - curr_logLik
  a3 <- log(pDEATH) - log(sum(unlist(lapply(prop_tree, '[[', 'NOG')))) - log(pBIRTH) + log(sum(unlist(lapply(curr_tree, '[[', 'l')))) 
  a4 <- gDEATH - gBIRTH
  a <- min(sum(a1, a2, a3, a4), 0)
  return(a)
}
rDEATH <- function(curr_tree, prop_tree, pBIRTH = 1/3, pDEATH = 1/3){
  
  # Identify pruned leaves and the proposed parent node (new leaf)
  curr_leaves <- which(lapply(curr_tree, length) > 0)
  prop_leaves <- which(lapply(prop_tree, length) > 0)
  pruned_leaves <- setdiff(curr_leaves, prop_leaves)
  bL <- pruned_leaves[1]; bR <- pruned_leaves[2]; b <- getParent(bL)
  
  # Retrieve the log-likelihood for the pruned leaves and the proposed parent node
  curr_logLik <- sum(vapply(c(bL, bR), function(l) curr_tree[[l]]$logLik, FUN.VALUE = 1))
  prop_logLik <- prop_tree[[b]]$logLik
  
  # Compute the splitting prior
  prop_rho <- pSplit(d = curr_tree[[b]]$d)
  curr_rho <- pSplit(d = curr_tree[[bL]]$d)
  
  # Compute the log-likelihood of the entire proposal
  gDEATH <- dnorm(prop_tree[[b]]$mu, prop_tree[[b]]$m, prop_tree[[b]]$v, log = TRUE)
  gBIRTH <- sum(unlist(lapply(curr_tree[c(bL, bR)], function(l) dnorm(l$mu, l$m, l$v, log = TRUE))))
  
  # Calculate the M-H acceptance ratio on the log scale
  a1 <- log(1 - prop_rho) - log(prop_rho) - 2 * log(1 - curr_rho)
  a2 <- prop_logLik - curr_logLik
  a3 <- log(pBIRTH) - log(sum(unlist(lapply(prop_tree, '[[', 'l')))) - log(pDEATH) + log(sum(unlist(lapply(curr_tree, '[[', 'NOG')))) # 1 if only root node in prop_tree 
  a4 <- gBIRTH - gDEATH
  a <- min(0, sum(a1, a2, a3, a4))
  
  return(a)
}
rCHANGE <- function(curr_tree, prop_tree){
  
  # Identify the leaf nodes affected by the change
  nodes <- which(lapply(curr_tree, length) > 0)
  bL <- nodes[(which(unlist(lapply(curr_tree, '[[', 'rule')) != unlist(lapply(prop_tree, '[[', 'rule'))) + 1)[1]]
  bR <- getSibling(bL)
  
  # Retrieve the log-likelihood under the current and proposed trees
  curr_logLik <- sum(vapply(c(bL, bR), function(l) curr_tree[[l]]$logLik, FUN.VALUE = 1))
  prop_logLik <- sum(vapply(c(bL, bR), function(l) prop_tree[[l]]$logLik, FUN.VALUE = 1))
  
  # Compute the log-likelihood of the proposal
  curr_gCHANGE <- sum(unlist(lapply(curr_tree[c(bL, bR)], function(l) dnorm(l$mu, l$m, l$v, log = TRUE))))
  prop_gCHANGE <- sum(unlist(lapply(prop_tree[c(bL, bR)], function(l) dnorm(l$mu, l$m, l$v, log = TRUE))))
  
  # Calculate the M-H acceptance ratio on the log scale
  a <- min(prop_logLik - curr_logLik + curr_gCHANGE - prop_gCHANGE, 0)
  return(a)
}
