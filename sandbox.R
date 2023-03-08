strata <- as.numeric(as.factor(strata))
windows <- table(strata)
system.time(table(strata))
system.time(tabulate(strata))
max.win <- max(windows)
na.locs <- max.win * unique(strata)[which(windows != max.win)]

n.obs <- length(unique(strata))

denoms <- NA * numeric(max.win*n.obs)

denoms[-na.locs] <- x %*% coef(m0)[1:3] + z * coef(m0)[4]
tmp_d <- matrix(denoms, ncol = max.win, byrow = TRUE)

logLik <- (denoms[-na.locs])[y == 1] - matrixStats::rowLogSumExps(tmp_d, na.rm = TRUE)


system.time(matrixStats::rowLogSumExps(tmp_d, na.rm = TRUE))

split(tmp_d, seq(1, ))

split(as.data.frame(tmp_d)[sample(1:nrow(tmp_d)),], rep(1:10, each = nrow(tmp_d)/10))

all(table(strata) == tabulate(strata))

system.time(table(as.numeric(as.factor(strata))))
system.time(tabulate(as.factor(strata)))


test <- as.matrix(data.frame(x1 = rep(1:2, each = 1e6), x2 = rep(1:2, times = 1e6)))
testmu <- as.matrix(data.frame(x1 = rep(1:2, each = 2), x2 = rep(1:2, times= 2), mu = 1:4))
test <- matrix(0, nrow = 1e4, ncol = 1e6)


mtest <- merge(test, testmu, by = colnames(test))
system.time(mtest <- merge(test, testmu, by = colnames(test))$mu)
system.time(mtest <- left_join(as.data.frame(test), as.data.frame(testmu), by = colnames(test))$mu)
test.dt <- data.table::as.data.table(test)
testmu.dt <- data.table::as.data.table(test)
data.table::setkey(test.dt, colnames(test))
system.time(mtest <- left_join(data.table::as.data.table(test), data.table::as.data.table(testmu), by = colnames(test)))


mu2 <- data.table::as.data.table(mu)
system.time(mu[2,] <- 1)
system.time(mu2[2,] <- 1)

