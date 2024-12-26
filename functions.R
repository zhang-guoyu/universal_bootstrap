library(ggplot2)
library(foreach)
library(parallel)
library(doParallel)
library(gridExtra)
library(reshape2)


pca <- function(k,X,Sigma,n,p) {
  A <- t(X) %*% (X)/n-Sigma
  result <- eigen(A)
  tk <- abs(result$values)[order(abs(result$values),decreasing = TRUE)[1:k]]
  max_eigenvalue <- sqrt(sum(tk*tk))
  max_eigenvalue
}

pcan <- function(X,Sigma,n,p) {
  Ta <- solve(msrqt(Sigma))
  Z <- X%*%Ta
  I <- diag(rep(1,p))
  A <- t(Z) %*% (Z)/n-I
  result <- eigen(A)
  max_eigenvalue <- max(abs(result$values))
  max_eigenvalue
}

pcann <- function(X,Sigma,n,p){
  sigma <- mean(diag(Sigma))
  t1 <- pca(1,X,Sigma,n,p)/sigma
  t2 <- pcan(X,Sigma,n,p)
  t1*t1+t2*t2
}




supn <- function(X,Sigma,n,p){
  #hat1x <- t(X) %*% (X)/n1
  #hat1y <- t(Y) %*% (Y)/n2
  #hat2x <- t(X*X) %*% (X*X)/n1-hat1x*hat1x
  #hat2y <- t(Y*Y) %*% (Y*Y)/n2-hat1y*hat1y
  Ta <- solve(msrqt(Sigma))
  Z <- X%*%Ta
  Z <- scale(Z, center=TRUE, scale=FALSE)
  I <- diag(rep(1,p))
  hat1 <- t(Z) %*% (Z)/n-I
  #hat2 <- t(Z*Z) %*% (Z*Z)/n-2*diag(colMeans(Z*Z))+I-hat1*hat1
  hat2 <- t(Z^2) %*% (Z^2) / n - (t(Z) %*% (Z)/n)^2
  R<-(hat1)/sqrt(hat2/n)
  #diag(R)<-0
  max(R^2)
}


flss <- function(X,Sigma,n,p){
  Ta <- solve(msrqt(Sigma))
  Z <- X%*%Ta
  Sn <- t(Z) %*% (Z)/n
  I <- diag(rep(1,p))
  mf <- mean(Z**4)
  W <- mean( diag( (Sn - I)%*%(Sn - I) ) )-p/n*(mean( diag(Sn))  )**2+p/n
  rfs <- n*W-p-(mf-2)
  rfs/2
}

utest <- function(X, Sigma,n,p){
  #n1 = nrow(X)
  #n2 = nrow(Y)
  #X.center = scale(X, center=TRUE, scale=FALSE)
  #Y.center = scale(Y, center=TRUE, scale=FALSE)
  Ta <- solve(msrqt(Sigma))
  Z <- X%*%Ta
  
  A = LC.Avalue(Z,n,p)
  C = LC.Cvalue(Z,n,p)
  Tn = A - 2*C + p
  #Tn.sd = A.x*2/n1 + A.y*2/n2
  
  #pVal = pnorm(Tn/Tn.sd, lower.tail=FALSE)
  #return(list(Tn = Tn, Tn.sd = Tn.sd, test.stat=Tn/Tn.sd, pVal = pVal))
  Tn
}

LC.Avalue <- function(X,n,p){
  n = nrow(X)
  sample.est <- (X %*% t(X))^2
  diag(sample.est) <- 0 ## exclude self-self products
  Avalue <- sum(sample.est) / (n * (n-1))
  
  return (Avalue)
}

LC.Cvalue <- function(X,n,p){
  #n1 = nrow(X)
  #n2 = nrow(Y)
  S <- t(X) %*% (X)/n
  Cvalue <- sum(diag(S))
  #sample.est <- (X %*% t(Y))^2 
  #Cvalue <- sum(sample.est) / (n1 * n2)
  
  return (Cvalue)
}


msrqt <- function(A){
  svd_result <- svd(A)
  mymatrix <- svd_result$u %*% sqrt(diag(svd_result$d)) %*% t(svd_result$v)
  mymatrix
}

mymatrixf <- function(method,decay,n,p){
  if (method == "AR1"){
    Sigma <- cbar *cbar *toeplitz(decay^(0:(p-1)))
    mymatrix <- msrqt(Sigma)
  }else if((method=="block")|(method=="block_nondiag")){
    block.num <- p/10
    off.value <- 0.55
    mymatrix <- matrix(0, nrow=p, ncol=p)
    block.size <- ceiling(p/block.num)
    for (k in c(1:block.num)) {
      i.start <- (k-1)*block.size + 1
      i.end <- min(k*block.size, p)
      mymatrix[i.start:i.end, i.start:i.end] <- off.value
    }
    diag(mymatrix) <- 1
    mymatrix <- msrqt(mymatrix)
  }else if((method=="eigendecay")|(method=="eigendecay_nondiag")   ){
    mymatrix <- diag(cbar*((1:(p))**(-decay/2)))
    mymatrix <- msrqt(mymatrix)
  }
  else if((method=="noise diag")   ){
    mymatrix <- diag(rexp(p))
    mymatrix <- msrqt(mymatrix)
  }
  else if((method=="inter AR1")   ){
    mymatrix <- matrix(nrow = p, ncol = p, data = sapply(1:p, function(i) {
      sapply(1:p, function(j) {
        (-1)^(i + j) * decay^abs(i - j)^0.5
      })
    }))
    mymatrix <- msrqt(mymatrix)
  }
  else{
    Sigma <- cbar *cbar *toeplitz(decay^(0:(p-1)))
    mymatrix <- msrqt(Sigma)
  }
  mymatrix
}


mymatrix2f <- function(mymatrix,dimethod,decay,sigma,n,p){
  Sigma <- mymatrix %*% mymatrix
  if (dimethod == "eigen"){
    svd_result <- svd(Sigma)
    myseq <- svd_result$d
    myseq[5] <- myseq[5]*sigma
    mymatrix2 <- svd_result$u %*% sqrt(diag(myseq)) %*% t(svd_result$v)
  }else if(dimethod=="nondiag"){
    re <- rnorm(p)
    mre <- sqrt(sum(re*re))
    re <- re/mre
    Sigma2 <- Sigma+sigma*re%*%t(re)
    
    re <- rnorm(p)
    mre <- sqrt(sum(re*re))
    re <- re/mre
    Sigma2 <- Sigma2+sigma*re%*%t(re)#matrix(sigma,p,p) 
    
    mymatrix2 <- msrqt(Sigma2)
  }else if(dimethod=="nondiag eigen"){
    svd_result <- svd(Sigma)
    myseq <- svd_result$d
    myseq[5] <- myseq[5]*(sigma/2+1)
    Sigma2 <- svd_result$u %*% diag(myseq) %*% t(svd_result$v)
    
    re <- rnorm(p)
    mre <- sqrt(sum(re*re))
    re <- re/mre
    Sigma2 <- Sigma2+sigma*re%*%t(re)/4#matrix(sigma,p,p) 
    
    mymatrix2 <- msrqt(Sigma2)
  }else if(dimethod=="block_diag"){
    block.num <- p
    off.value <- sigma
    Sigma2 <- matrix(0, nrow=p, ncol=p)
    block.size <- ceiling(p/block.num)
    for (k in c(1:block.num)) {
      i.start <- (k-1)*block.size + 1
      i.end <- min(k*block.size, p)
      Sigma2[i.start:i.end, i.start:i.end] <- off.value*0.45
    }
    diag(Sigma2) <- off.value
    Sigma2 <- Sigma2+Sigma
    mymatrix2 <- msrqt(Sigma2)
  }
  else{
    mymatrix2 <- mymatrix
  }
  mymatrix2
}