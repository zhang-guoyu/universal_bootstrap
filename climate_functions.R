#####################################################################
# fingerprint via tls                                               #
#####################################################################
## X: the fingerprinting patterns
## Y: observed responses
## nruns.X: number of runs for the fingerprinting
## cov: covariance matrix estimate from Z1
## Z.2: separate group of control runs used in variance estimate
## conf.level: confidence level
## TS: indicating whether using the two sample approach or the proposed method
tlsFingerprint <- function(X, Y, nruns.X, cov, Z.2, precision = FALSE, conf.level = 0.90, TS = TRUE) {
  ## Z.2 is the second sample for the variance estimation
  X <- as.matrix(X)
  if (! precision) {
    tmpMat <- eigen(cov)
    cov.sinv <- Re(tmpMat$vectors %*% diag(1 / sqrt(tmpMat$values)) %*% 
                     t(tmpMat$vectors))
    ## cov.sinv <- Re(sqrtm(cov))
  } else {
    tmpMat <- eigen(cov)
    cov.sinv <- Re(tmpMat$vectors %*% diag(sqrt(tmpMat$values)) %*% 
                     t(tmpMat$vectors))
  }
  
  ## compute the estimation using internal function tlsLm
  output <- tlsLmTS(cov.sinv %*% X, cov.sinv %*% Y, nruns.X, conf.level, Z.2, cov.sinv, TS = TS)
  beta.hat <- output$beta.hat
  ci.estim <- output$ci
  sd.estim <- output$sd

  ## label the cols and rows of the estimation
  rownames(ci.estim) <- c("B ANT", "B NAT",
                          "N ANT", "N NAT")
  colnames(ci.estim) <- c("lower", "upper")
  rownames(sd.estim) <- c("B", "N")
  colnames(sd.estim) <- c("ANT", "NAT")
  
  names(beta.hat) <- c("ANT", "NAT")
  ## combine the results
  result <- list(coefficient = beta.hat,
                 confidence.interval = ci.estim, 
                 var.est = sd.estim)
  result
}

## estimate via Total least square approach
tlsLmTS <- function(X, Y, nruns.X, conf.level, Z.2, cov.sinv, TS) {
  ## input: 
  ##   X: n*k matrix, including k predictors
  ##   Y: n*1 matrix, the observations
  ##   nruns.X: the number of runs used for computing each columns of X
  ##   Z.2: the second sample to estimate the variance
  ##   cov.sinv: weight matrix
  ## output: 
  ##   a list containing the estimate and confidence interval of the scaling 
  ##   factors estimate.
  ##   coefficient: a k vector, the scaling factors best-estimates
  ##   confidence interval: the lower and upper bounds of the confidence 
  ##   interval on each scaling factor. 
  ##   dcons: the variable used in the residual consistency check. 
  ##   X.tilde: a k*n matrix, the reconstructed responses patterns TLS fit,
  ##   Y.tilde: a 1*n matrix, the reconstructed observations TLS fit.
  if (! is.matrix(Y)) {
    stop("Y should be a n*1 matrix")
  }
  if (dim(X)[1] != dim(Y)[1])  {  ## check size of X and Y
    stop("sizes of inputs X, Y are not consistent")
  }
  n <- dim(Y)[1]  ## number of observations
  m <- dim(X)[2]  ## number of predictors
  ## Normalise the variance of X
  X <- X * t(sqrt(nruns.X) %*% matrix(1, 1, n))
  if(length(nruns.X) != 1) {
    Dn.X <- diag(sqrt(nruns.X))
  } else {
    Dn.X <- sqrt(nruns.X)
  }
  #### for the proposed one sample approach
  Estls <- function(X, Y, Dn.X) {
    M <- cbind(X, Y)
    
    lambda <- tail(eigen(t(M) %*% M)$values, 1)
    
  }
  
  #### for the two sample approach
  Estls <- function(X, Y, Dn.X, TS) {
    M <- cbind(X, Y)
    ## get the smallest eigenvalue, used for variance estimate
    lambda <- tail(eigen(t(M) %*% M)$values, 1)
    ## beta.hat for the tls regression for adjusted X and Y with equal variance
    beta.hat1 <- as.vector(solve(t(X) %*% X - lambda * diag(m)) %*% t(X) %*% Y)
    
    if(TS) {
      ## if use two step approach with separate group of control runs Z2
      #### update the eigenvalues and corresponding covariance estimate based
      #### on Ribes et al. 2013
      Lambda2 <- diag(t(svd(M)$u) %*% M %*% t(M) %*% svd(M)$u) / 
        diag(t(svd(M)$u) %*% cov.sinv %*% cov(Z.2) %*% cov.sinv %*% svd(M)$u)
      XtX <- (svd(M)$v %*% diag(Lambda2) %*% t(svd(M)$v))[1:m, 1:m]
      Delta.hat2 <- (XtX - tail(Lambda2, 1) * diag(m)) / n
      sigma2.hat2 <- tail(Lambda2, 1) / n
      ## [I|beta.hat1]
      I.b <- cbind(diag(m), beta.hat1)
      ## var.hat for beta.hat1
      Var.hat1 <- sigma2.hat2 * (1 + sum(beta.hat1^2)) *
        solve(Delta.hat2) %*% (Delta.hat2 + sigma2.hat2 *  solve(I.b %*% t(I.b))) %*% solve(Delta.hat2) / n
    } else {
      sigma2.hat <- lambda / n
      Delta.hat <- (t(X) %*% X - lambda * diag(m)) / n
      ## beta.hat1 <- as.vector(RegGTLS(Y, X))
      ## [I|beta.hat1]
      I.b <- cbind(diag(m), beta.hat1)
      ## var.hat for beta.hat1
      Var.hat1 <- sigma2.hat * (1 + sum(beta.hat1^2)) * 
        (solve(Delta.hat) + sigma2.hat * solve(Delta.hat) %*% 
           solve(I.b %*% t(I.b)) %*% solve(Delta.hat)) / n
    }
    
    ## beta.hat and var.hat for the un prewhitening X and Y
    beta.hat <- beta.hat1 %*% Dn.X
    Var.hat <- diag(Var.hat1) * nruns.X
    ## return the results
    list(beta.hat = beta.hat, Var.hat = Var.hat)
  }
  
  #### function for the bootstrap
  Estls.beta <- function(X, Y, Dn.X) {
    M <- cbind(X, Y)
    lambda <- tail(eigen(t(M) %*% M)$values, 1)
    ## beta.hat for the tls regression for adjusted X and Y with equal variance
    beta.hat1 <- as.vector(solve(t(X) %*% X - lambda * diag(m)) %*% t(X) %*% Y)
    ## beta.hat and var.hat for the un prewhitening X and Y
    beta.hat1 %*% Dn.X
  }
  
  ## compute the estimation
  tmp.res <- Estls(X, Y, Dn.X, TS = TS)
  beta.hat <- tmp.res$beta.hat
  var.hat <- tmp.res$Var.hat
  
  ## get the normal critical value
  Z.crt <- qnorm((1 - conf.level) / 2, lower.tail = FALSE)
  ## confidence interval from asymptotical normal distribution
  ci.norm <- cbind(t(beta.hat - Z.crt * sqrt(var.hat)), 
                   t(beta.hat + Z.crt * sqrt(var.hat)))
  colnames(ci.norm) <- c(0.5 - conf.level / 2, 0.5 + conf.level / 2)
  
  ## compute the normal standard deviation
  sd.norm <- sqrt(var.hat)
  
  ## compute the confidence interval from bootstrap
  B <- 1000
  ## nonparamatric bootstrap
  ## error control with tryCatch function
  for(i in 1:10) {
    beta.s <- tryCatch({
      resample <- sapply(1:B,
                         function(x) {
                           sample(1:n, size = n, replace = TRUE)
                         })
      beta.s <- apply(resample, 2,
                      function(x) {
                        Xs <- X[x, ]
                        Ys <- Y[x, ]
                        Estls.beta(Xs, Ys, Dn.X)
                      })
      if(ncol(X) == 1) {
        matrix(beta.s, ncol = B)
      } else {
        beta.s
      }
    }, error = function(e) {
      as.matrix(c(0, 0))
    })
    if (ncol(beta.s) == B) {
      break
    }
  }
  
  ## alpha value
  alpha <- 1 - conf.level
  
  ## compute the bootstrap confidence interval
  ci.ordboot <- t(apply(beta.s, 1, 
                        function(x) {
                          x <- rmout(x)
                          quantile(x, c(alpha / 2, alpha / 2 + conf.level))
                  }))
  ## compute the bootstrap standard deviation
  sd.ordboot <- t(apply(beta.s, 1, 
                      function(x) {
                        x <- rmout(x)
                        sd(x)
                      }))
  
  list(beta.hat = beta.hat, ci = rbind(ci.ordboot, ci.norm), 
       sd = rbind(sd.ordboot, sd.norm))
}

##############################################################################
## tuning covariance                                                        ##
##############################################################################
## function for calibration simulation procedure
## parameters:
##   X: response to the external forcings
##   Y: real world observations 
##   nruns.X: numer of runs for each external forcing
##   cov: estimated covariance matrix (assumed to be the substitution of truth in calibration simulation)
##   Z.2: the second sample to estimate the covariance matrix
##   n.sample: sample size for the control runs
##   B: number of replicates for the simulation
tuningCov <- function(X, Y, nruns.X, cov, n.sample, B = 500, TS = FALSE) {
  X <- as.matrix(X)
  ## solve the inverse sqrt of covariance matrix
  tmpMat <- eigen(cov)
  cov.sinv <- Re(tmpMat$vectors %*% diag(1 / sqrt(tmpMat$values)) %*% 
                     t(tmpMat$vectors))
  ## internal function to get the expected X and beta for simulation
  tls <- function(X, Y, nruns.X) {
    if (! is.matrix(Y)) {
      stop("Y should be a n*1 matrix")
    }
    if (dim(X)[1] != dim(Y)[1])  {  ## check size of X and Y
      stop("sizes of inputs X, Y are not consistent")
    }
    n <- dim(Y)[1]  ## number of observations
    m <- dim(X)[2]  ## number of predictors
    ## Normalise the variance of X
    X <- X * t(sqrt(nruns.X) %*% matrix(1, 1, n))
    Estls <- function(X, Y) {
      M <- cbind(X, Y)
      svd.M <- svd(M)
      output <- svd.M$u %*% diag(c(svd.M$d[1:m], 0)) %*% t(svd.M$v)
      colnames(output) <- c("ANT", "NAT", "Y")
      output
    }
    output <- Estls(X, Y)
    output[, c("ANT", "NAT")] <- output[, c("ANT", "NAT")] / 
      t(sqrt(nruns.X) %*% matrix(1, 1, n))
    output
  }
  ## compute the expected X and beta for simulation
  output <- tls(cov.sinv %*% X, cov.sinv %*% Y, nruns.X)
  output <- Re(tmpMat$vectors %*% diag(sqrt(tmpMat$values)) %*% 
               t(tmpMat$vectors)) %*% output
  ## calibration simulation procedure
  ## save the list of estimated scaling factors beta and corresponding variances  
  out.beta <- out.var <- NULL
  for (i in 1:B) {
    ## use tryCatch and loop to avoid error in improper data generation
    for(er in 1:10) {
      tmp.new <- tryCatch({
        ## generate new bootstrap data set based on the estimated expected responses
        Y.new <- MASS::mvrnorm(n = 1, mu = output[, "Y"], Sigma = cov)
        X.new <- cbind(MASS::mvrnorm(n = 1, mu = output[, "ANT"], Sigma = cov / nruns.X[1]), 
                       MASS::mvrnorm(n = 1, mu = output[, "NAT"], Sigma = cov / nruns.X[2]))
        ## generate the set of control runs
        ctlruns.new <- genClt(n.sample, 1, cov)
        Cov.new <- Covest(ctlruns.new[[1]][1:(n.sample/2), ], loss = "l2")$output
        ## use the fingerprinting function to refit the model and get estimations of beta and variances
        tlsFingerprint(X.new, Y.new, nruns.X, Z.2 = ctlruns.new[[1]][-(1:(n.sample/2)), ], Cov.new, TS = TS)
      }, error = function(e) {
        ""
      })
      if (length(tmp.new) > 1) {
        break
      }
    }
    if(length(tmp.new) > 1) {
      out.beta <- rbind(out.beta, tmp.new$coefficient)
      out.var <- rbind(out.var, as.vector(tmp.new$var.est))
    }
  }
  list(rbind(rowMeans(apply((out.beta), 2, sd) / t(rmoutInt(out.var[, c(1, 3)]))), 
             rowMeans(apply((out.beta), 2, sd) / t(rmoutInt(out.var[, c(2, 4)])))), 
       out.beta, out.var)
  # list(out.beta, out.var)
}