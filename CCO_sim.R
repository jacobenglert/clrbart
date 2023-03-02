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

# Exposure (zero-mean wave)
set.seed(144)
z <- rnorm(length(t), sin(2*pi*as.numeric(t)/3) + cos(2*pi*as.numeric(t)/3), sd = 1) # rnorm(length(t), 0, 1)
# plot(z ~ t); lines(z ~ t)
z <- ifelse(z > 2, 1, 0) # convert exposure to binary
exposure <- data.frame(Date = t, Z = z)
Z <- rep(z, each = n)

# Subject-level covariates
x_probs <- rep(0.5, 5) # population frequencies
set.seed(8451)
x <- vapply(x_probs, function(p) rbinom(n, 1, p), FUN.VALUE = 1:n)
X <- matrix(rep(t(x), length(t)), ncol = ncol(x), byrow = TRUE)
Xc <- 1 - X
colnames(Xc) <- paste0('X', 1:ncol(X),'c')

# Set up design matrix
# XZ <- cbind(Z, X[,1]*Z, X[,2]*Z, X[,1]*X[,2]*Z)
XZ <- Z * cbind(Xc[,1]*Xc[,2], Xc[,1]*X[,2], X[,1]*Xc[,2], X[,1]*X[,2])

# Simulate cases based on pre-specified odds ratios
beta <- matrix(c(0.7, 0.7, -0.4, 0.1), ncol = 1)
alpha <- rep(-7, n) #rnorm(n, mean = -10, sd = 1)
expit <- function(x){exp(x)/(1+exp(x))}
p.OR <- expit(rep(alpha, each = length(t)) + XZ %*% beta)
# p.RR <- exp(rep(alpha, each = length(t)) + XZ %*% beta)
set.seed(2187)
y <- rbinom(length(p.OR), 1, p.OR)

# Identify exposures and covariates associated with the cases only
y.star <- y[y == 1]
t.star <- rep(t, each = n)[y == 1]
X.star <- matrix(X[y == 1,], ncol = ncol(x), byrow = FALSE)
colnames(X.star) <- paste0('X', 1:ncol(X.star))

# Create CC structure
CC <- make_cc(t, z, t.star, X.star)

# Recover the odds ratios
# m0 <- clogit(Y ~ Z + X1:Z + X2:Z + X1:X2:Z + strata(ID), data = CC)
# summary(m0)
# logLik(m0)

# Recover exact odds ratios for subsets
Xc <- matrix(!select(CC, starts_with('X')), ncol = ncol(X)) * 1
colnames(Xc) <- paste0(colnames(X.star),'c')
CC2 <- cbind(CC, Xc)
m01 <- clogit(Y ~ X1c:X2c:Z + X1c:X2:Z + X1:X2c:Z + X1:X2:Z + strata(ID), data = CC2)
m02 <- clogit(Y ~ X1c:Z + X1:X2c:Z + X1:X2:Z + strata(ID), data = CC2)
summary(m02)
logLik(m02)

library(clogitLasso)


write_csv(CC, 'Data/test.csv')

XCL <- as.matrix(as.vector(CC2$Z) * as.matrix(select(CC2, starts_with('X'))))
YCL <- as.vector(CC2$Y)
IDCL <- as.vector(CC2$ID)
m03 <- clogitLasso(XCL, YCL, IDCL)
