# Program Name: CCO_sim.R
# Author: Jacob Englert
# Date: 11 January 2023
# Purpose: Simulate conditional logistic regression data with multiple binary
# predictors.


# Load Packages -----------------------------------------------------------
library(tidyverse)
library(lubridate)
library(survival)
source('make_cc.R')

# Simulation --------------------------------------------------------------
# Number of subjects
n <- 100000

# Time points (days)
t <- seq(ymd('2022-01-01'), ymd('2022-12-31'), by = 'days')

# Time-varying Exposure (zero-mean wave)
set.seed(144)
z <- rnorm(length(t), sin(2*pi*as.numeric(t)/3) + cos(2*pi*as.numeric(t)/3), sd = 1) # rnorm(length(t), 0, 1)
# plot(z ~ t); lines(z ~ t)
z <- ifelse(z > 2, 1, 0) # convert exposure to binary
Z <- rep(z, each = n)

# Time-varying confounders
set.seed(854)
x <- replicate(3, rnorm(length(t)))
colnames(x) <- paste0('X', 1:ncol(x))
X <- x %x% rep(1, n) # do.call(rbind, replicate(n, x, simplify = FALSE))

# Subject-level effect modifiers
w_probs <- rep(0.5, 5) # population frequencies
set.seed(8451)
w <- vapply(w_probs, function(p) rbinom(n, 1, p), FUN.VALUE = 1:n)
colnames(w) <- paste0('W', 1:ncol(w))
W <- do.call(rbind, replicate(length(t), w, simplify = FALSE)) # w %x% rep(1, length(t))
Wc <- 1 - W
colnames(Wc) <- paste0('W', 1:ncol(W),'c')


# Set up design matrix
dim(Z); dim(X); dim(W)
XWZ <- cbind(X,  Z * cbind(Wc[,1]*Wc[,2], Wc[,1]*W[,2], W[,1]*Wc[,2], W[,1]*W[,2]))

# Simulate cases based on pre-specified odds ratios
beta <- matrix(c(-0.5, 0.5, 0.1, 0.7, 0.7, -0.4, 0.1), ncol = 1)
alpha <- rep(-8, n) #rnorm(n, mean = -10, sd = 1)
expit <- function(x){exp(x)/(1+exp(x))}
p.OR <- expit(rep(alpha, each = length(t)) + XWZ %*% beta) 
# p.RR <- exp(rep(alpha, each = length(t)) + XZ %*% beta)
set.seed(2187)
Y <- rbinom(length(p.OR), 1, p.OR)

# Identify exposures, confounders, and effect modifiers associated with the cases only
Y.star <- Y[Y == 1]
t.star <- rep(t, each = n)[Y == 1]
W.star <- matrix(W[Y == 1,], ncol = ncol(W), byrow = FALSE)
colnames(W.star) <- paste0('W', 1:ncol(W.star))

# Create CC structure
CC <- make_cc(t, cbind(x, Z = z), t.star, W.star)

# Recover exact odds ratios for subsets
Wc.star <- matrix(!select(CC, starts_with('W')), ncol = ncol(W)) * 1
colnames(Wc.star) <- paste0(colnames(W.star),'c')
CC2 <- cbind(CC, Wc.star)
mWZ <- clogit(Y ~ W1c:W2c:Z + W1c:W2:Z + W1:W2c:Z + W1:W2:Z + strata(ID), data = CC2)
mXWZ <- update(mWZ, ~ . + X1 + X2 + X3)
summary(mWZ)
summary(mXWZ)

mWZ2 <- clogit(Y ~ W1c:Z + W1:W2c:Z + W1:W2:Z + strata(ID), data = CC2)
mXWZ2 <- update(mWZ2, ~ . + X1 + X2 + X3)
summary(mWZ2)
summary(mXWZ2)

logLik(mWZ)
logLik(mXWZ)
logLik(mWZ2)
logLik(mXWZ2)

# Export test dataset
write_csv(CC, 'Data/XWZ.csv')

