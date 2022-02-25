rm(list=ls())

start=Sys.time()
format(start,"%H:%M:%S")

alpha.t <- numeric() 

T2 <- numeric()
T_1 <- numeric() ; T_2.1 <- numeric() ; T_3.1_2 <- numeric()


## original rho
rho123 <- c(-0.3,0,0.3,0.6)

## corresponding rho
rho12 <- c(-0.40,0,0.40,0.80) ; 
rho13 <- c(-0.31,0,0.31,0.62) ; 
rho23 <- c(-0.42,0,0.42,0.78) 


num <- c(20,30,40,50)


## in-control parameters
mu <- 0 ; sig <- 1
p0 <- 0.2
lam0 <- 3


## Load pacakge
library(readr)
library(mvtnorm)

d_df_final <- NULL

for(d in 1:length(rho123)){
  
  c_df_final <- NULL
  
  for(k in 1:length(num)){
    
    N <- 100000
    for(j in 1:N){
  
      ## Setting m (observations) ; n (subgroup size)
      n <- num[k]  
  
      CorrZ <- matrix(0,3,3) # CorrZ will be the correlation matrix of the standard normal vector.
      CorrZ[1,2] <- rho12[d]      
      CorrZ[1,3] <- rho13[d]
      CorrZ[2,3] <- rho23[d]
      CorrZ <- CorrZ + t(CorrZ) 
      diag(CorrZ)<-1
  
      muZ <- c(0,0,0)   
      
      library(mvtnorm)
      Z <- rmvnorm(n,muZ,CorrZ)
      
      U <- pnorm(Z)
      
      Y <- matrix( cbind(qnorm(U[,1],mu,sig),qbinom(U[,2],1,p0),qpois(U[,3],lam0)),n,3)
      
      Y_bar <- matrix(colMeans(Y),3,1)
      
      CorrY <- matrix(0,3,3) 
      CorrY[1,2] <- rho123[d]      
      CorrY[1,3] <- rho123[d]
      CorrY[2,3] <- rho123[d]
      CorrY <- CorrY + t(CorrY)          
      diag(CorrY) <- c(sig,p0*(1-p0),lam0)
      
      mu_Y <- matrix(c(mu,p0,lam0),3,1) ; S_Y <- matrix(CorrY,3,3)
      
      ## Compute T^2
      invS <- solve(S_Y)
      T2[j] <- n*t(Y_bar-mu_Y)%*% invS %*% (Y_bar-mu_Y)
      
      
      ## Decomposition start
      
      ## Compute T_1
      T_1[j] <- n*(Y_bar[1,1]-mu_Y[1,1])^2 / S_Y[1,1]
      
      ## Compute T_2.1
      b2 <- S_Y[1,2] / S_Y[1,1] 
      X_2.1 <- mu_Y[2,1] + b2 * (Y_bar[1,1]-mu_Y[1,1])
      s_2.1 <- S_Y[2,2] - (S_Y[1,2])*(S_Y[1,1])^(-1)*(S_Y[1,2])
      T_2.1[j] <- n*(Y_bar[2,1]-X_2.1)^2 / s_2.1
      
      
      ## Compute T_3.1,2
      Sxx <- matrix(c(S_Y[1,1],S_Y[1,2],S_Y[2,1],S_Y[2,2]),2,2)
      Sx <- matrix(c(S_Y[1,3],S_Y[2,3]),2,1) 
      b3 <- solve(Sxx) %*% Sx 
      x_2 <- matrix(c(Y_bar[1,1],Y_bar[2,1]),2,1) ; x_2_bar <- matrix(c(mu_Y[1,1],mu_Y[2,1]),2,1)
      X_3.1 <- mu_Y[3,1] + t(b3) %*% (x_2 - x_2_bar )
      s_3.1 <- S_Y[3,3] - t(Sx) %*% solve(Sxx) %*% (Sx)
      T_3.1_2[j] <- n*(Y_bar[3,1]-X_3.1)^2 / s_3.1
      
      
      ##--- T2 chart ---------------
      alpha.t[j] <- round(1 - pchisq(T2[j],3),5)

    }  
    
    alpha <- 0.005

    ##==========  Type I ===================================
    Type1_B <- sum(alpha.t < alpha)/N
    #Type1_B 
    
    c_df <- cbind(num[k],Type1_B)
    c_df <- data.frame(c_df)
    names(c_df)[1] <- c("n") 
    
    c_df_temp <- c_df
    c_df_final  <- rbind(c_df_final,c_df_temp)
    
  }
  
  rho <- matrix(rho1234[d],k,1)
  rho <- data.frame(rho)
  names(rho) <- c("rho")
  
  d_df <- cbind(rho,c_df_final)
  d_df <- data.frame(d_df)
  
  d_df_temp <- d_df
  d_df_final  <- rbind(d_df_final,d_df_temp) 


}


d_df_final

write_csv(d_df_final,file=paste("Result/p03/Hotelling-T2-Type1.csv",sep = ""),
          append = TRUE,col_names = TRUE)



end=Sys.time()
end-start