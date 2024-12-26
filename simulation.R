library(ggplot2)
library(foreach)
library(parallel)
library(doParallel)
library(gridExtra)
library(reshape2)



set.seed(100)

setwd("")

#source('sLED.R', encoding = 'UTF-8')
#source('sLEDpermute.R', encoding = 'UTF-8')
source('symmPMD.R', encoding = 'UTF-8')
source("functions.R")

sz=Sys.time()
print(sz)

cl <- makeCluster(100)
registerDoParallel(cl)


cbar <- 1
nlist <- 300
plist <- 1000
size <- 0.95

dimethod <- "nondiag eigen" #c("nondiag eigen" ) #"nondiag eigen"# "nondiag"     "eigen"  "nondiag eigen" "block_diag" "null"
methodl <-   c("block")#c("AR1","block","noise diag") # "eigendecay""block""noise diag" "AR1" "inter AR1"
distril <- c("Gau")#c("t","unif","Gau m t","Gau m unif","t m unif")




num <- 2000#number of repitation
bnum <- 3000#number of bootstrap




lambda_seq <- seq(0,3,0.3)#(seq(0,2,0.2))#sqrt(seq(1,10,0.5))#sqrt(seq(0,5,0.2))   #sqrt(seq(1,1.2,0.01)) for AR1 heavy tails, sqrt(seq(1,30,2)) for AR1 elliptic and elliptic
#sqrt(seq(1,60,3)) for  heavy tails, sqrt(seq(1,10,0.5)) for AR1
#sqrt(seq(0,10,0.5))/p for AR1 nondiag
#sqrt(seq(1,5,0.2)) for block, sqrt(seq(0,5,0.2)) for nondiag block
trial <- length(lambda_seq)

decay_seq <- 0.6#c(0,0.4,0.6)#c(0.4,0.6,1,1.2)#c(0,0.2,seq(0.4,0.8,0.1))
case <- length(decay_seq)



prob_pca1 <- rnorm(trial)
prob_pcan <- rnorm(trial)
prob_pcann <- rnorm(trial)
prob_supn <- rnorm(trial)
prob_flss <- rnorm(trial)
prob_utest <- rnorm(trial)



for(distri in distril){
  for(n in nlist){
    for(method in methodl){
      for(p in plist){
        for (alpha in c(1:case)){
          for (m in c(1:trial)){
            buff <-paste("_one samples_bootstrap_test_table_N_",n,"_p_",p,"_",method,"_",dimethod,"_",distri)
            
            
            
            print(buff)
            s=Sys.time()
            print(s)
            
            cat("p:", p, "\n")
            cat("n:", n, "\n")
            decay <- decay_seq[alpha]
            cat("decay:", decay, "\n")
            sigma <- lambda_seq[m]
            cat("lambda:", sigma, "\n")
            
            cat("method:", method, "\n")
            cat("distri:", distri, "\n")
            
            
            mymatrix <- mymatrixf(method = method,decay = decay,n,p)
            mymatrix2 <- mymatrix2f(mymatrix = mymatrix,dimethod = dimethod,decay = decay,sigma=sigma,n,p)
            Sigma <- mymatrix%*%mymatrix
            
           
            
            
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
            
            threshhold_pca1 <- quantile(thresh_pca1,size)
            threshhold_pca2 <- quantile(thresh_pca2,size)
            threshhold_pca5 <- quantile(thresh_pca5,size)
            threshhold_pcan <- quantile(thresh_pcan,size)
            threshhold_pcann <- quantile(thresh_pcann,size)
            threshhold_supn <- quantile(thresh_supn,size)
            
            
            testresult <- foreach( i = c(1:num),.combine = rbind,.errorhandling = "remove")%dopar%{
              #print(i)
              set.seed(i+1000)
              
              if(distri == "nonGau"){
                a1 <- rnorm(n*p)
                a2 <- rnorm(n*p)
                b1 <- seq(from=0.5,to=1,length.out=n*p)
                b2 <- seq(from=1,to=2,length.out=n*p)
                x <- (b1*(a1*a1-1)+b2*(a2*a2-1))/sqrt(2*(b1*b1+b2*b2))
                Y <- matrix(x , nrow = n)
              }else if(distri == "t"){
                df <- 12
                Y <- matrix(rt(n*p,df=df) , nrow = n)/sqrt(df/(df-2))
              }else if(distri == "unif"){
                Y <- matrix(runif(n*p,-1,1) , nrow = n)*sqrt(3)
              }else if(distri == "Gau m t"){
                Y1 <- matrix(rnorm(n*p/2) , nrow = n)
                df <- 12
                Y2 <- matrix(rt(n*p/2,df=df) , nrow = n)/sqrt(df/(df-2))
                Y <- cbind(Y1, Y2)
              }else if(distri == "Gau m unif"){
                Y1 <- matrix(rnorm(n*p/2) , nrow = n)
                Y2 <- matrix(runif(n*p/2,-1,1) , nrow = n)*sqrt(3)
                Y <- cbind(Y1, Y2)
              }else if(distri == "t m unif"){
                df <- 12
                Y1 <- matrix(rt(n*p/2,df=df) , nrow = n)/sqrt(df/(df-2))
                Y2 <- matrix(runif(n*p/2,-1,1) , nrow = n)*sqrt(3)
                Y <- cbind(Y1, Y2)
              }else{
                Y <- matrix(rnorm(n*p) , nrow = n)
              }
              X <- Y%*%mymatrix2
              
              
              
              test_pca1 <- pca(1,X,Sigma,n,p)
              threshhold_pcaub1 <- threshhold_pca1
              testre_pca1 <- (test_pca1>threshhold_pcaub1)
              
        
              
              test_pcan <- pcan(X,Sigma,n,p)
              threshhold_pcan <- threshhold_pcan
              testre_pcan <- (test_pcan>threshhold_pcan)
              
              test_pcann <- pcann(X,Sigma,n,p)
              threshhold_pcann <- threshhold_pcann
              testre_pcann <- (test_pcann>threshhold_pcann)
              
              test_supn <- supn(X,Sigma,n,p)
              threshhold_supn <- (4 * log(p) - log(log(p)) - log(8 * pi) - 2 * log(log((size)^(-1))))#supnb(X,Sigma)#threshhold_supn##supnb(X,Y)
              testre_supn <- (test_supn>threshhold_supn)
              
              test_flss <- flss(X,Sigma,n,p)
              threshhold_flss <- qnorm(size)
              testre_flss <- (test_flss>threshhold_flss)
              
              test_utest <- utest(X,Sigma,n,p)
              threshhold_utest <- 2*qnorm(size)*sqrt( p*(p+1) / (n*(n-1)) )#utestb(X,Sigma)
              testre_utest <- (test_utest>threshhold_utest)
              
              c(test_pca1,test_pca2,test_pca5,test_pcan,test_pcann,test_supn,test_flss,test_utest,
                testre_pca1,testre_pca2,testre_pca5,testre_pcan,testre_pcann,testre_supn,testre_flss,testre_utest)
            }
            
            
            
            test_pca1 <- testresult[,1]
            test_pca2 <- testresult[,2]
            test_pca5 <- testresult[,3]
            test_pcan <- testresult[,4]
            test_pcann <- testresult[,5]
            test_supn <- testresult[,6]
            test_flss <- testresult[,7]
            test_utest <- testresult[,8]
            
            testre_pca1 <- testresult[,9]
            testre_pca2 <- testresult[,10]
            testre_pca5 <- testresult[,11]
            testre_pcan <- testresult[,12]
            testre_pcann <- testresult[,13]
            testre_supn <- testresult[,14]
            testre_flss <- testresult[,15]
            testre_utest <- testresult[,16]
            
            
            
            prob_pca1[m] <- sum(testre_pca1)/num
            prob_pca2[m] <- sum(testre_pca2)/num
            prob_pca5[m] <- sum(testre_pca5)/num
            prob_pcan[m] <- sum(testre_pcan)/num
            prob_pcann[m] <- sum(testre_pcann)/num
            prob_supn[m] <- sum(testre_supn)/num
            prob_flss[m] <- sum(testre_flss)/num
            prob_utest[m] <- sum(testre_utest)/num
            
            cat("reject probability of supn:", prob_supn[m], "\n")
            cat("reject probability of flss:", prob_flss[m], "\n")
            cat("reject probability of utest:", prob_utest[m], "\n")
            cat("reject probability of pca1:", prob_pca1[m], "\n")
            cat("reject probability of pca2:", prob_pca2[m], "\n")
            cat("reject probability of pca5:", prob_pca5[m], "\n")
            cat("reject probability of pcan:", prob_pcan[m], "\n")
            cat("reject probability of pcann:", prob_pcann[m], "\n")
            
            
            e=Sys.time()
            print(e-s)
          }
          write.table(cbind(lambda_seq,prob_pca1), file = paste("decay_",decay,"/bootstrap",buff,"_pca1_power_curve.txt"), sep = ",",quote = FALSE)
          write.table(cbind(lambda_seq,prob_pca2), file = paste("decay_",decay,"/bootstrap",buff,"_pca2_power_curve.txt"), sep = ",",quote = FALSE)
          write.table(cbind(lambda_seq,prob_pca5), file = paste("decay_",decay,"/bootstrap",buff,"_pca5_power_curve.txt"), sep = ",",quote = FALSE)
          write.table(cbind(lambda_seq,prob_pcan), file = paste("decay_",decay,"/bootstrap",buff,"_pcan_power_curve.txt"), sep = ",",quote = FALSE)
          write.table(cbind(lambda_seq,prob_pcann), file = paste("decay_",decay,"/bootstrap",buff,"_pcann_power_curve.txt"), sep = ",",quote = FALSE)
          write.table(cbind(lambda_seq,prob_supn), file = paste("decay_",decay,"/bootstrap",buff,"_supn_power_curve.txt"), sep = ",",quote = FALSE)
          write.table(cbind(lambda_seq,prob_flss), file = paste("decay_",decay,"/bootstrap",buff,"_flss_power_curve.txt"), sep = ",",quote = FALSE)
          write.table(cbind(lambda_seq,prob_utest), file = paste("decay_",decay,"/bootstrap",buff,"_utest_power_curve.txt"), sep = ",",quote = FALSE)
          
          probf <- data.frame(
            lambda_seq,
            prob_supn,
            prob_flss,
            prob_utest,
            prob_pcan,
            prob_pca1,
            prob_pcann
          )
          
          probmelted <- melt(probf,id="lambda_seq")
          colnames(probmelted)[2] = 'Method'
          plot_1 <- ggplot(data = probmelted,aes(x=lambda_seq,y=value,group = Method,
                                                 color=Method,shape=Method))+
            geom_line(size=1.2)+
            geom_point(size=3)+
            scale_shape_manual(values = c(1:6))+
            geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
            theme_bw() +
            scale_shape_manual(values = c(1:6), labels = c('Supn','Lfn','Ufn','Roy','Opn','Com')) +
            scale_colour_discrete(labels=c('Supn','Lfn','Ufn','Roy','Opn','Com'))+
            scale_x_continuous(name="Signal Level")+
            scale_y_continuous(name="Empirical Power")+
            theme(panel.grid.minor = element_blank(),legend.position = c(.08,.83),
                  legend.box.background = element_rect(color="black"),
                  axis.title.x=element_text(size=44),
                  axis.title.y=element_text(size=44),
                  axis.text.x=element_text(size=40),
                  axis.text.y=element_text(size=40),
                  legend.title =element_text(size=28),
                  legend.text = element_text(size=28)) 
          ggsave(paste("decay_",decay,"/bootstrap",buff,"_power_curve.png"),plot=plot_1,width=16,height=12)
          
          
          
        }
      }
    }
  }
}






stopImplicitCluster()
stopCluster(cl)

ez=Sys.time()
print(ez-sz)

