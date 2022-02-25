rm(list=ls())

start=Sys.time()
format(start,"%H:%M:%S")

alpha.x <- numeric() ; alpha.p <- numeric() ; alpha.c <- numeric()

zbar <- numeric() ; zp <- numeric() ; zcbar <- numeric()

## original rho
rho123 <- c(-0.3,0,0.3,0.6)

## corresponding rho
rho12 <- c(-0.40,0,0.40,0.80) ; 
rho13 <- c(-0.31,0,0.31,0.62) ; 
rho23 <- c(-0.42,0,0.42,0.78) 

num <- c(20,30,40,50)

## Load pacakge
library(readr)
library(mvtnorm)

d_df_final <- NULL

for(d in 1:length(rho123)){
  
  c_df_final <- NULL
  
  for(k in 1:length(num)){
    
    N <- 500000
    for(j in 1:N){
      
      ## Setting n (observations) 
      n <- num[k]  
      
      CorrZ <- matrix(0,3,3) # CorrZ will be the correlation matrix of the standard normal vector.
      CorrZ[1,2] <- rho12[d]      
      CorrZ[1,3] <- rho13[d]
      CorrZ[2,3] <- rho23[d]
      CorrZ <- CorrZ + t(CorrZ)          # Get the lower triangle too.
      diag(CorrZ)<-1
      
      muZ <- c(0,0,0)   
      Z <- rmvnorm(n,muZ,CorrZ)
      
      U <- pnorm(Z) 
      
      Y <- matrix( cbind(qnorm(U[,1]), qbinom(U[,2],1,0.3) , qpois(U[,3],3)), n,3   )

      
      ##--- X-bar chart ---------
      mu <- 0 ; sig <- 1
      
      z <- Y[,1]
      zbar[j] <- mean(z) 
      
      lower.x <- pnorm((zbar[j]-mu)/(sig/sqrt(n)))
      upper.x <- ( 1- pnorm((zbar[j]-mu)/(sig/sqrt(n))) )
      
      alpha.x[j] <- round(2*min(c(lower.x,upper.x)),5)
      
   
      ##--- p-chart ---------
      p0 <- 0.3
      
      zp[j] <- sum(Y[,2])/n
      
      lower.p <- pnorm((zp[j]-p0)/sqrt(p0*(1-p0)/n))
      upper.p <- ( 1- pnorm( (zp[j]-p0)/sqrt(p0*(1-p0)/n) ) )
      
      alpha.p[j] <- round(2*min(c(lower.p,upper.p)),5)
      
   
      ##--- c-chart ---------
      lambda <- 3
      
      zc <- Y[,3]
      zcbar[j] <- mean(zc)
      
      lower.c <- pnorm((zcbar[j]-lambda)/sqrt(lambda/n))
      upper.c <- ( 1- pnorm((zcbar[j]-lambda)/sqrt(lambda/n)) )
      
      alpha.c[j] <- round(2*min(c(lower.c,upper.c)),5)
      
      
    }  
    
    alpha <- 0.005

    ##==========  Type I ===================================
    A <- as.matrix(cbind(alpha.x,alpha.p,alpha.c),N,3)
    
    ## Bonferroni's 
    a <- matrix(0,N,1) 
    A_ord <- matrix(0,N,3) 
    for(i in 1:N){
      A_ord[i,] <- sort(A[i,])
      if(A_ord[i,1] < alpha/3 || A_ord[i,2] < alpha/3 || A_ord[i,3] < alpha/3 ){
        a[i] <- 1 
      }
    }
    Type1_B <- sum(a)/N
    #Type1_B
    
    
    ## Holm step-down procedure
    A_ord <- matrix(0,N,3) 
    aH <- matrix(0,N,1) 
    for(i in 1:N){
      A_ord[i,] <- sort(A[i,])
      if(A_ord[i,1] < alpha/3 || A_ord[i,2] < alpha/2 || A_ord[i,3] < alpha){
        aH[i] <- 1 
      }
    }
    Type1_H <- sum(aH)/N
    #Type1_H
    
    
    c_df <- cbind(num[k],Type1_B,Type1_H)
    c_df <- data.frame(c_df)
    names(c_df) <- c("n","Type1_B","Type1_H") 
    
    c_df_temp <- c_df
    c_df_final  <- rbind(c_df_final,c_df_temp)
    
  }
  
  rho <- matrix(rho123[d],k,1)
  rho <- data.frame(rho)
  names(rho) <- c("rho")
  
  d_df <- cbind(rho,c_df_final)
  d_df <- data.frame(d_df)
  
  d_df_temp <- d_df
  d_df_final  <- rbind(d_df_final,d_df_temp) 


}


d_df_final

write_csv(d_df_final,file=paste("Result/p03/Holm-and-Bon-Type1.csv",sep = ""),
          append = TRUE,col_names = TRUE)



end=Sys.time()
end-start