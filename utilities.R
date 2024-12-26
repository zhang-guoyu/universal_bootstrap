##############################################################################
# Covariance estimate via Linear shrinkage estimator                         #
##############################################################################

## covariance matrix estimator under frobenius (l2) norm
Covest <- function(X, loss = c("l2", "stein", "mv")) {
  ## the other two estimators are not included in this manuscript
  loss <- match.arg(loss)
  ## check the method
  if (!loss %in% c("l2", "stein", "mv")) {
    stop("please select one objective function from 'l2', 'stein' or 'mv'")
  }
  output <- NULL
  if (loss == "l2") {
    ## Linear shrinkage estimator under frobenious loss function
    output <- lwRegcov(X)
  } else {
    ## not included in this manuscript
    stop("type of estimator is not included")
  }
  list(output = output)
}

## compute the regularized estimate of the covariance matrix
## References:
## Ledoit O., M. Wolf (2004) A well-conditioned estimator for
## large-dimensional covariance matrices.
## Journal of Multivariate Analysis, 88(2): 365-411.
## L&W estimate

lwRegcov <- function(X) {
  ## input:
  ##   the independent replicates of the random vector, denoted by matrix X
  ##   each column represents a random variable. Thus X is n*p matrix
  ## output:
  ##   regularized covariance matrix
  n <- nrow(X)
  p <- ncol(X)
  sample.cov <- cov(X)
  Ip <- diag(p)
  m <- sum(diag(sample.cov %*% Ip)) / p ## first estimate in L&W
  Xp <- sample.cov - m * Ip
  d2 <- sum(diag(Xp %*% t(Xp))) / p  ## second estimate in L&W
  
  bt <- (diag(X %*% t(X))^2 - 2 * diag(X %*% sample.cov %*% t(X)) + rep(1, n) *
    sum(diag(sample.cov %*% sample.cov))) / p
  
  bb2 <- 1 / n^2 * sum(bt)
  b2 <- min(bb2,d2)  ## third estimate in L&W
  a2 <- d2 - b2  ## fourth estimate in L&W
  ## the regularized estimate of the covariance matrix
  ## b2 * m / d2 * Ip + a2 / d2 * sample.cov
  eigen(sample.cov)$vector %*% 
  diag(b2 * m / d2 + a2 / d2 * eigen(sample.cov)$values) %*% 
  t(eigen(sample.cov)$vector)
}

##############################################################################
# internal functions                                                         #
##############################################################################

## functions for removing possible outliers
rmout <- function(x) {
  ind <- c(which(x %in% boxplot.stats(x)$out),  which(is.na(x)))
  ind <- unique(ind)
  if(length(ind) == 0) {
    x
  } else {
    x[-ind]
  }
}

rmoutInt <- function(x) {
  ind <- c(which(x[, 1] %in% boxplot.stats(x[, 1])$out), 
           which(x[, 2] %in% boxplot.stats(x[, 2])$out), 
           which(is.na(x[, 1])), 
           which(is.na(x[, 2])))
  ind <- unique(ind)
  if(length(ind) == 0) {
    x
  } else {
    x[-ind, , drop = FALSE]
  }
}

## function for generating cltruns for given sample size
## n: sample size
## B: number of replicates
## Cov: covariance matrix
genClt <- function(n, B, Cov) {
  lapply(1:B, 
         function(x) {
           MASS::mvrnorm(n, mu = rep(0, nrow(Cov)), Sigma = Cov)
         })
}