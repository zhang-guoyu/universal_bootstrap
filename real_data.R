library(ggplot2)
library(foreach)
library(parallel)
library(doParallel)
library(gridExtra)
library(reshape2)

library(methods)

setwd("")
set.seed(419)

source("functions.R")
source("utilities.R")
source("climate_functions.R")

load("appDat.rdata")

## set the output dir
outDir <- "res/"
#### get the name of the regions
list.file <- names(dat)

name <- "global_40x30(GL).rdata"
cat(name, "\n")

size <- 0.95
bnum <- 10000




s=Sys.time()
print(s)

cl <- makeCluster(100)
registerDoParallel(cl)


for (ite in c(1:5)){
  cat('year',1950+ite*10,'-',1960+ite*10, "\n")
  ytr <- 54*2
  len <- length(dat[[name]]$observation)
  veck <- c((len-(6-ite)*ytr) : (len-(5-ite)*ytr))
  observation <- dat[[name]]$observation[veck]#[1:ytr]
  ANT.forcing <- dat[[name]]$ANT.forcing[veck,]#[1:ytr]
  NAT.forcing <- dat[[name]]$NAT.forcing[veck,]#[1:ytr]
  ctlruns <- dat[[name]]$ctlruns[veck,]#[1:ytr]
  
  ind <- which(! is.na(observation))
  
  ## remove missing values
  ctlruns <- ctlruns[ind, ]
  ANT.forcing <- ANT.forcing[ind, ]
  NAT.forcing <- NAT.forcing[ind, ]
  
  X1 <- scale(t(ANT.forcing), center=TRUE, scale=FALSE)
  X2 <- scale(t(NAT.forcing), center=TRUE, scale=FALSE)
  X = rbind(X1,X2)
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  eps <- scale(t(ctlruns), center=TRUE, scale=FALSE)
  kn <- 0
  Xe <- eps[1:kn,]
  
  noi <- 2.5
  X <- X + noi*matrix(rnorm(n*p) , nrow = n)
  X <- scale(X, center=TRUE, scale=FALSE)
  #eps <- eps + noi*matrix(rnorm(dim(eps)[1]*dim(eps)[2]) , nrow = dim(eps)[1])
  #eps <- scale(eps, center=TRUE, scale=FALSE)
  
  
  Sigma <- t(eps[(kn+1):dim(eps)[1],])%*%eps[(kn+1):dim(eps)[1],]/(dim(eps)[1]-kn)
  Sigma <- Sigma + diag(rep(noi*2,p))
  
  
  
  
  cat("p:", p, "\n")
  cat("n:", n, "\n")
  mymatrix <- msrqt(Sigma)
  
  
  threshresult <- foreach( i = c(1:bnum),.combine = rbind,.errorhandling = "remove")%dopar%{
    #print(i)
    set.seed(i+1000)
    Z <- matrix(rnorm(n*p) , nrow = n)
    X <- Z%*%mymatrix
    I <- diag(rep(1,p))
    
    
    A <- t(X) %*% (X)/n-Sigma
    B <- t(Z) %*% (Z)/n-I
    
    resulta <- eigen(A)
    resultb <- eigen(B)
    
    t1 <- abs(resulta$values)[order(abs(resulta$values),decreasing = TRUE)[1]]
    thresh_pca1 <- sqrt(sum(t1*t1))
    
    t2 <- abs(resulta$values)[order(abs(resulta$values),decreasing = TRUE)[1:2]]
    thresh_pca2 <- sqrt(sum(t2*t2))
    
    t5 <- abs(resulta$values)[order(abs(resulta$values),decreasing = TRUE)[1:5]]
    thresh_pca5 <- sqrt(sum(t5*t5))
    
    tn1 <- max(abs(resultb$values))
    thresh_pcan <- tn1
    
    sigma <- mean(diag(Sigma))
    thresh_pcann <- (t1/sigma)^2 +tn1^2
    
    hat <- t(Z*Z) %*% (Z*Z)/n-2*diag(colMeans(Z*Z))+I-B*B
    R<-(B)/(sqrt(hat/n))
    thresh_supn <- max(R*R)
    
    c(thresh_pca1,thresh_pca2,thresh_pca5,thresh_pcan,thresh_pcann,thresh_supn,
      thresh_pca1,thresh_pca2,thresh_pca5,thresh_pcan,thresh_pcann,thresh_supn)
  }
  
  thresh_pca1 <- threshresult[,1]
  thresh_pca2 <- threshresult[,2]
  thresh_pca5 <- threshresult[,3]
  thresh_pcan <- threshresult[,4]
  thresh_pcann <- threshresult[,5]
  thresh_supn <- threshresult[,6]
  
  test_pca1 <- pca(1,X,Sigma,n,p)
  testre_pca1 <- 1-mean(test_pca1>thresh_pca1)
  
  test_pca2 <- pca(2,X,Sigma,n,p)
  testre_pca2 <- 1-mean(test_pca2>thresh_pca2)
  
  test_pca5 <- pca(5,X,Sigma,n,p)
  testre_pca5 <- 1-mean(test_pca5>thresh_pca5)
  
  test_pcan <- pcan(X,Sigma,n,p)
  testre_pcan <- 1-mean(test_pcan>thresh_pcan)
  
  test_pcann <- pcann(X,Sigma,n,p)
  testre_pcann <- 1-mean(test_pcann>thresh_pcann)
  
  test_supn <- supn(X,Sigma,n,p)
  testre_supn <-  1-exp(-exp( (4 * log(p) - log(log(p)) - log(8 * pi) - test_supn )/2   )) 
  
  test_flss <- flss(X,Sigma,n,p)
  testre_flss <- 1-pnorm(test_flss)
  
  test_utest <- utest(X,Sigma,n,p)
  testre_utest <-  (1-pnorm( (test_utest)/(2*sqrt( p*(p+1) / (n*(n-1)) ))    ) )  
  
  cat('treat group','\n')
  cat("p value of supn:", testre_supn, "\n")
  cat("p value of flss:", testre_flss, "\n")
  cat("p value of utest:", testre_utest, "\n")
  cat("p value of pca1:", testre_pca1, "\n")
  #cat("p value of pca2:", testre_pca2, "\n")
  cat("p value of pcan:", testre_pcan, "\n")
  cat("p value of pcann:", testre_pcann, "\n")
  
  
  
  Z <- matrix(rnorm(n*p) , nrow = n)
  Xe <- Z%*%mymatrix
  
  
  test_pca1e <- pca(1,Xe,Sigma,n,p)
  testre_pca1 <- 1-mean(test_pca1e>thresh_pca1)
  
  test_pca2e <- pca(2,Xe,Sigma,n,p)
  testre_pca2 <- 1-mean(test_pca2e>thresh_pca2)
  
  test_pca5e <- pca(5,Xe,Sigma,n,p)
  testre_pca5 <- 1-mean(test_pca5e>thresh_pca5)
  
  test_pcane <- pcan(Xe,Sigma,n,p)
  testre_pcan <- 1-mean(test_pcane>thresh_pcan)
  
  test_pcanne <- pcann(Xe,Sigma,n,p)
  testre_pcann <- 1-mean(test_pcanne>thresh_pcann)
  
  test_supne <- supn(Xe,Sigma,n,p)
  testre_supn <-  1-exp(-exp( (4 * log(p) - log(log(p)) - log(8 * pi) - test_supne )/2   )) 
  
  test_flsse <- flss(Xe,Sigma,n,p)
  testre_flss <- 1-pnorm(test_flsse)
  
  test_uteste <- utest(Xe,Sigma,n,p)
  testre_utest <-  (1-pnorm( (test_uteste)/(2*sqrt( p*(p+1) / (n*(n-1)) ))    ) )  
  
  cat('control group','\n')
  cat("p value of supn:", testre_supn, "\n")
  cat("p value of flss:", testre_flss, "\n")
  cat("p value of utest:", testre_utest, "\n")
  cat("p value of pca1:", testre_pca1, "\n")
  #cat("p value of pca2:", testre_pca2, "\n")
  cat("p value of pcan:", testre_pcan, "\n")
  cat("p value of pcann:", testre_pcann, "\n")
  ez=Sys.time()
  print(ez-s)
}



stopImplicitCluster()
stopCluster(cl)

ez=Sys.time()
print(ez-s)


