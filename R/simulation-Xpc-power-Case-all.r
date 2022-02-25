rm(list=ls())

start=Sys.time()
format(start,"%H:%M:%S")

alpha.x <- numeric() ; alpha.p <- numeric() ; alpha.c <- numeric()
pval.x <- numeric() ; pval.p <- numeric() ; pval.c <- numeric()

zbar <- numeric() ; zp <- numeric() ; zcbar <- numeric()
zbar1 <- numeric() ; zp1 <- numeric() ; zcbar1 <- numeric()

T2 <- numeric()
T_1 <- numeric() ; T_2.1 <- numeric() ; T_3.1_2 <- numeric()

T12 <- numeric()
T1_1 <- numeric() ; T1_2.1 <- numeric() ; T1_3.1_2 <- numeric()

alpha.t <- numeric() ; pval.t <- numeric()

alpha.Tx <- numeric() ; alpha.Tp <- numeric() ; alpha.Tc <- numeric()
pval.Tx <- numeric() ; pval.Tp <- numeric() ; pval.Tc <- numeric()

## in-control parameters
mu <- 0 ; sig <- 1
p0 <- 0.3
lambda <- 3

## out-of-control parameters
s <- c(0.2,0.3,0.4,0.5)

sn <- sig*s
mu1_v <- round(c(mu,mu+sn),2)  

sp <- c(0.05,0.10,0.15,0.20)
p1_v <- round(c(p0,p0+sp),2)   

sc <- c(0.2,0.4,0.6,0.8)
lambda1_v <- round(c(lambda,lambda+sc),2)  

## original rho
rho123 <- c(-0.3,0,0.3,0.6)
rho1234 <- c(-0.3,0,0.3,0.6)

## corresponding rho
rho12 <- c(-0.40,0,0.40,0.80) ; 
rho13 <- c(-0.31,0,0.31,0.62) ; 
rho23 <- c(-0.42,0,0.42,0.78) 

## Holm and Bom Adjust alpha
a20 <- c(0.0058,0.0059,0.0054,0.0059)
a30 <- c(0.0044,0.0044,0.0048,0.0048)
a40 <- c(0.0057,0.0057,0.0057,0.0057)
a50 <- c(0.0058,0.0058,0.0058,0.0058)

aBH <- matrix(rbind(a20,a40,a30,a50),4,4)

## Hotelling Adjust alpha
aT20  <- c(0.0001,0.0053,0.0002,0.0606)  # n20
aT30  <- c(0.0001,0.0054,0.0002,0.0596)  # n30
aT40  <- c(0.0001,0.0050,0.0002,0.0618)  # n40
aT50  <- c(0.0001,0.0054,0.0002,0.0601)  # n50


aT <- matrix(rbind(aT20,aT40,aT30,aT50),4,4)
num <- c(20,40,30,50)

## Load pacakge
library(readr)
library(mvtnorm)

d_df_final <- NULL

for(nu in 1:length(num)){
  
  for(d in 1:4){
    c_df_final <- NULL 
  
    for(s in 1:5){
      for(l in 1:5){
        for (k in 1:5) {
        
        start_time <- proc.time()

        N <- 100000
        for(j in 1:N){
          
          ## Setting n (observations) 
          n <- num[nu]  
          
          CorrZ <- matrix(0,3,3) # CorrZ will be the correlation matrix of the standard normal vector.
          CorrZ[1,2] <- rho12[d]      
          CorrZ[1,3] <- rho13[d]
          CorrZ[2,3] <- rho23[d]
          CorrZ <- CorrZ + t(CorrZ)          # Get the lower triangle too.
          diag(CorrZ)<-1
          
          muZ <- c(0,0,0)   
          Z <- rmvnorm(n,muZ,CorrZ)
          
          U <- pnorm(Z) 
          
          Y <- matrix( cbind(qnorm(U[,1],mu,sig), qbinom(U[,2],1,p0) , qpois(U[,3],lambda)), n,3   )
          
          mu1 <- mu1_v[s] ; p1 <- p1_v[l] ; lambda1 <- lambda1_v[k]
          Y1 <- matrix( cbind(qnorm(U[,1],mu1,sig), qbinom(U[,2],1,p1) , qpois(U[,3],lambda1)), n,3   )
          
          
          ##--- X-bar chart ---------
          #mu <- 0 ; sig <- 1
          
          z <- Y[,1]
          zbar[j] <- mean(z) 
          
          lower.x <- pnorm((zbar[j]-mu)/(sig/sqrt(n)))
          upper.x <- ( 1- pnorm((zbar[j]-mu)/(sig/sqrt(n))) )
          
          alpha.x[j] <- round(2*min(c(lower.x,upper.x)),5)
          
          ## H1 : mu1
          z1 <- Y1[,1]
          zbar1 <- mean(z1)
          
          lower1.x <- pnorm((zbar1-mu)/(sig/sqrt(n)))
          upper1.x <- ( 1- pnorm((zbar1-mu)/(sig/sqrt(n))) )
          
          pval.x[j] <- round(2*min(c(lower1.x,upper1.x)),5)
          
          ##--- p-chart ---------
          #p0 <- 0.3
          
          zp[j] <- sum(Y[,2])/n
          
          lower.p <- pnorm((zp[j]-p0)/sqrt(p0*(1-p0)/n))
          upper.p <- ( 1- pnorm( (zp[j]-p0)/sqrt(p0*(1-p0)/n) ) )
          
          alpha.p[j] <- round(2*min(c(lower.p,upper.p)),5)
          
       
          ## H1 : p1
          zp1[j] <- sum(Y1[,2])/n
          
          lower1.p <- pnorm((zp1[j]-p0)/sqrt(p0*(1-p0)/n))
          upper1.p <- ( 1- pnorm((zp1[j]-p0)/sqrt(p0*(1-p0)/n)) )
          
          pval.p[j] <- round(2*min(c(lower1.p,upper1.p)),5)
          
          ##--- c-chart ---------
          #lambda <- 3
          
          zc <- Y[,3]
          zcbar[j] <- mean(zc)
          
          lower.c <- pnorm((zcbar[j]-lambda)/sqrt(lambda/n))
          upper.c <- ( 1- pnorm((zcbar[j]-lambda)/sqrt(lambda/n)) )
          
          alpha.c[j] <- round(2*min(c(lower.c,upper.c)),5)
          
          ## H1 : lambda1
          zc1 <- Y1[,3]
          zcbar1[j] <- mean(zc1)
          
          lower1.c <- pnorm((zcbar1[j]-lambda)/sqrt(lambda/n))
          upper1.c <- ( 1- pnorm((zcbar1[j]-lambda)/sqrt(lambda/n)) )
          
          pval.c[j] <- round(2*min(c(lower1.c,upper1.c)),5)  
          
          
          
          ##=============== Hotelling T2 chart =================================
          Y_bar <- matrix(colMeans(Y),3,1)
          
          CorrY <- matrix(0,3,3) 
          CorrY[1,2] <- rho123[d]      
          CorrY[1,3] <- rho123[d]
          CorrY[2,3] <- rho123[d]
          CorrY <- CorrY + t(CorrY)          
          diag(CorrY) <- c(sig,p0*(1-p0),lambda)
          
          mu_Y <- matrix(c(mu,p0,lambda),3,1) ; S_Y <- matrix(CorrY,3,3)
          
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
          
          ##--- X-bar chart ---------
          alpha.Tx[j] <- round(1 - pchisq(T_1[j],1),5)
        
          ##--- p-chart ---------
          alpha.Tp[j] <- round(1 - pchisq(T_2.1[j],1),5)    
          
          ##--- c-chart ---------
          alpha.Tc[j] <- round(1 - pchisq(T_3.1_2[j],1),5)
          
          
          ## H1
          Y1_bar <- matrix(colMeans(Y1),3,1)
          
          ## Compute T^2
          invS <- solve(S_Y)
          T12[j] <- n*t(Y1_bar-mu_Y)%*% invS %*% (Y1_bar-mu_Y)
          
          
          ## Decomposition start
          
          ## Compute T_1
          T1_1[j] <- n*(Y1_bar[1,1]-mu_Y[1,1])^2 / S_Y[1,1]
          
          ## Compute T_2.1
          b2 <- S_Y[1,2] / S_Y[1,1] 
          X_2.1 <- mu_Y[2,1] + b2 * (Y1_bar[1,1]-mu_Y[1,1])
          s_2.1 <- S_Y[2,2] - (S_Y[1,2])*(S_Y[1,1])^(-1)*(S_Y[1,2])
          T1_2.1[j] <- n*(Y1_bar[2,1]-X_2.1)^2 / s_2.1
          
          
          ## Compute T_3.1,2
          Sxx <- matrix(c(S_Y[1,1],S_Y[1,2],S_Y[2,1],S_Y[2,2]),2,2)
          Sx <- matrix(c(S_Y[1,3],S_Y[2,3]),2,1) 
          b3 <- solve(Sxx) %*% Sx 
          x_2 <- matrix(c(Y1_bar[1,1],Y1_bar[2,1]),2,1) ; x_2_bar <- matrix(c(mu_Y[1,1],mu_Y[2,1]),2,1)
          X_3.1 <- mu_Y[3,1] + t(b3) %*% (x_2 - x_2_bar )
          s_3.1 <- S_Y[3,3] - t(Sx) %*% solve(Sxx) %*% (Sx)
          T1_3.1_2[j] <- n*(Y1_bar[3,1]-X_3.1)^2 / s_3.1
          
          
          ##--- T2 chart ---------------
          pval.t[j] <- round(1 - pchisq(T12[j],3),5)
          
          ## mu1
          pval.Tx[j] <- round(1 - pchisq(T1_1[j],1),5)
          
          ##  p1
          pval.Tp[j] <- round(1 - pchisq(T1_2.1[j],1),5)
          
          ##  lambda1
          pval.Tc[j] <- round(1 - pchisq(T1_3.1_2[j],1),5) 
          
        }  
        
        alpha <- aBH[nu,d] ; alpha_T <- aT[nu,d]
        
        ############################################################
        ##======================  Type I ===================================
        A <- as.matrix(cbind(alpha.x,alpha.p,alpha.c),N,3)
        
        B1 <- matrix(0,N,3) 
        for(i in 1:N){
          for(j in 1:3){
            if(A[i,j] < alpha/3){ B1[i,j] <- 1 }
          }
        }
        
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
        Type1_B
        
        
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
        Type1_H
        
        ## Hotelling T2 chart
        Type1_TB <- sum(alpha.t < alpha_T)/N
        Type1_TB 
        
        ############################################################
        ##======================= Type II =========================
        
        P <- as.matrix(cbind(pval.x,pval.p,pval.c),N,3)
        
        ##--------- Bonferroni's Test ------------------------------------
        p <- matrix(0,N,1) 
        P_ord <- matrix(0,N,3) 
        for(i in 1:N){
          P_ord[i,] <- sort(P[i,])
          if(P_ord[i,1] < alpha/3 || P_ord[i,2] < alpha/3 || P_ord[i,3] < alpha/3 ){
            p[i] <- 1 
          }
        }
        Power_B <- sum(p)/N
        Power_B
        
        ## If any of three p-value < alpha/3, set to 1
        B11 <- matrix(0,N,3) 
        for(i in 1:N){
          for(j in 1:3){
            if(P[i,j] < alpha/3){ B11[i,j] <- 1 }
          }
        }
        
        p.mat <- B11
        
        c.ind_1 <- which(rowSums(p.mat)==1)
        c.ind_2 <- which(rowSums(p.mat)==2)
        c.ind_3 <- which(rowSums(p.mat)==3)
        
        c.ind_1_p1 <- length(which(p.mat[c.ind_1,1]==1))
        c.ind_1_p2 <- length(which(p.mat[c.ind_1,2]==1))
        c.ind_1_p3 <- length(which(p.mat[c.ind_1,3]==1))
        
        c.ind_2_p12 <- length(which(p.mat[c.ind_2,1]==1 & p.mat[c.ind_2,2]==1))
        c.ind_2_p13 <- length(which(p.mat[c.ind_2,1]==1 & p.mat[c.ind_2,3]==1))
        c.ind_2_p23 <- length(which(p.mat[c.ind_2,2]==1 & p.mat[c.ind_2,3]==1))
        
        c.ind_3_p123 <- length(which(p.mat[c.ind_3,1]==1 & p.mat[c.ind_3,2]==1 &
                                       p.mat[c.ind_3,3]==1))
    
        
        ##--------- Holm step-down procedure ----------------------------------
        P_ord <- matrix(0,N,3) 
        pH <- matrix(0,N,1)
        for(i in 1:N){
          P_ord[i,] <- sort(P[i,])
          if(P_ord[i,1] < alpha/3 || P_ord[i,2] < alpha/2 || P_ord[i,3] < alpha){
            pH[i] <- 1 
          }
        }
        
        Power_H <- sum(pH)/N
        Power_H
        
        P_ord <- matrix(0,N,3) ; P_sort <- matrix(0,N,3) 
        pH_indx <- matrix(0,N,3)
        for(i in 1:N){
          ## position in sort
          P_ord[i,] <- order(P[i,])
          ## sort
          P_sort[i,] <- as.matrix(sort(P[i,]),1,3)
          
          if(P_sort[i,1] < alpha/3){
            pH_indx[i,P_ord[i,1]] <- 1
          }
          if(P_sort[i,2] < alpha/2){
            pH_indx[i,P_ord[i,2]] <- 1
          }
          if(P_sort[i,3] < alpha){
            pH_indx[i,P_ord[i,3]] <- 1
          }
        }
        
        c.ind_H_1 <- which(rowSums(pH_indx)==1)
        c.ind_H_2 <- which(rowSums(pH_indx)==2)
        c.ind_H_3 <- which(rowSums(pH_indx)==3)
        
        c.ind_H_1_p1 <- length(which(pH_indx[c.ind_H_1,1]==1))
        c.ind_H_1_p2 <- length(which(pH_indx[c.ind_H_1,2]==1))
        c.ind_H_1_p3 <- length(which(pH_indx[c.ind_H_1,3]==1))
        
        c.ind_H_2_p12 <- length(which(pH_indx[c.ind_H_2,1]==1 & pH_indx[c.ind_H_2,2]==1))
        c.ind_H_2_p13 <- length(which(pH_indx[c.ind_H_2,1]==1 & pH_indx[c.ind_H_2,3]==1))
        c.ind_H_2_p23 <- length(which(pH_indx[c.ind_H_2,2]==1 & pH_indx[c.ind_H_2,3]==1))
        
        c.ind_H_3_p123 <- length(which(pH_indx[c.ind_H_3,1]==1 & pH_indx[c.ind_H_3,2]==1 
                                       & pH_indx[c.ind_H_3,3]==1))
        
        
        ##--------- Hotelling T2 chart ---------------------------------------
        ind_T <- which(pval.t < alpha_T)
        Power_TB <- sum(pval.t < alpha_T)/N
        Power_TB 
        
        TP <- as.matrix(cbind(pval.Tx,pval.Tp,pval.Tc),N,3)
        TP <- data.frame(TP)
        if(length(ind_T) != 0){
          TP_res <- TP[ind_T,]
        }else{
          TP_res <- matrix(0,1,3)
        }
        
        ## 
        B_T <- matrix(0,length(TP_res[,1]),3) 
        for(i in 1:length(TP_res[,1])){
          for(j in 1:3){
            if(TP_res[i,j] < alpha_T/3){ B_T[i,j] <- 1 }
          }
        }
        
        Tp.mat <- B_T
        
        c.ind_T_1 <- which(rowSums(Tp.mat)==1)
        c.ind_T_2 <- which(rowSums(Tp.mat)==2)
        c.ind_T_3 <- which(rowSums(Tp.mat)==3)
        
        c.ind_T_1_p1 <- length(which(Tp.mat[c.ind_T_1,1]==1))
        c.ind_T_1_p2 <- length(which(Tp.mat[c.ind_T_1,2]==1))
        c.ind_T_1_p3 <- length(which(Tp.mat[c.ind_T_1,3]==1))
        
        c.ind_T_2_p12 <- length(which(Tp.mat[c.ind_T_2,1]==1 & Tp.mat[c.ind_T_2,2]==1))
        c.ind_T_2_p13 <- length(which(Tp.mat[c.ind_T_2,1]==1 & Tp.mat[c.ind_T_2,3]==1))
        c.ind_T_2_p23 <- length(which(Tp.mat[c.ind_T_2,2]==1 & Tp.mat[c.ind_T_2,3]==1))
        
        c.ind_T_3_p123 <- length(which(Tp.mat[c.ind_T_3,1]==1 & Tp.mat[c.ind_T_3,2]==1 
                                     & Tp.mat[c.ind_T_3,3]==1))
        
        
        c_df_B <- c(c.ind_1_p1,c.ind_1_p2,c.ind_1_p3,
                    c.ind_2_p12,c.ind_2_p13,c.ind_2_p23,c.ind_3_p123)
        
        c_df_H <- c(c.ind_H_1_p1,c.ind_H_1_p2,c.ind_H_1_p3,
                    c.ind_H_2_p12,c.ind_H_2_p13,c.ind_H_2_p23,c.ind_H_3_p123)
        
        c_df_T <- c(c.ind_T_1_p1,c.ind_T_1_p2,c.ind_T_1_p3,
                    c.ind_T_2_p12,c.ind_T_2_p13,c.ind_T_2_p23,c.ind_T_3_p123)
        
        ## Calculate probability of correct diagnostics
        x1 <- s-1 ; bp1 <- l-1 ; c1 <- k-1
        pcd_B <- 0 ; pcd_H <- 0 ; pcd_T <- 0
        if(x1 != 0 & bp1 == 0 & c1 == 0){   ## Case 1 : only Normal
          pcd_B <- c.ind_1_p1/sum(p)
          pcd_H <- c.ind_H_1_p1/sum(pH)
          pcd_T <- c.ind_T_1_p1/length(ind_T)
        }
        if(x1 == 0 & bp1 != 0 & c1 == 0){  ## Case 2 : only Binomial
          pcd_B <- c.ind_1_p2/sum(p)
          pcd_H <- c.ind_H_1_p2/sum(pH)
          pcd_T <- c.ind_T_1_p2/length(ind_T)
        }
        if(x1 == 0 & bp1 == 0 & c1 != 0){  ## Case 3 : only Poisson
          pcd_B <- c.ind_1_p3/sum(p)
          pcd_H <- c.ind_H_1_p3/sum(pH)
          pcd_T <- c.ind_T_1_p3/length(ind_T)
        }
        
        if(x1 != 0 & bp1 != 0 & c1 == 0){   ## Case 4 :  Normal and Binomial
          pcd_B <- c.ind_2_p12/sum(p)
          pcd_H <- c.ind_H_2_p12/sum(pH)
          pcd_T <- c.ind_T_2_p12/length(ind_T)
        }
        if(x1 != 0 & bp1 == 0 & c1 != 0){  ## Case 5 : Normal and Poisson
          pcd_B <- c.ind_2_p13/sum(p)
          pcd_H <- c.ind_H_2_p13/sum(pH)
          pcd_T <- c.ind_T_2_p13/length(ind_T)
        }
        if(x1 == 0 & bp1 != 0 & c1 != 0){  ## Case 6 : Binomial and Poisson
          pcd_B <- c.ind_2_p23/sum(p)
          pcd_H <- c.ind_H_2_p23/sum(pH)
          pcd_T <- c.ind_T_2_p23/length(ind_T)
        }
        
        if(x1 != 0 & bp1 != 0 & c1 != 0){  ## Case 7 : all
          pcd_B <- c.ind_3_p123/sum(p)
          pcd_H <- c.ind_H_3_p123/sum(pH)
          pcd_T <- c.ind_T_3_p123/length(ind_T)
        }
        
        if(pcd_T == "NaN"){ pcd_T <- 0 }
        
        time <- proc.time() - start_time
        time_elapse <- time[3]
        
        c_df <- matrix(c(Type1_H,Type1_B,Type1_TB,Power_H,Power_B,Power_TB,
                         pcd_H,pcd_B,pcd_T,time_elapse),1,10)
        c_df <- data.frame(c_df)
        names(c_df) <- c("Type1_H","Type1_B","Type1_TB","Power_H","Power_B","Power_TB",
                         "pcd_H","pcd_B","pcd_T","time_elapse")
        
        
        count <- matrix(c(rho1234[d],x1,bp1,c1),1,4)
        count <- data.frame(count)
        names(count) <- c("rho","x","p","c")
        
        c_df <- cbind(count,c_df)
        
      
        ## save result to csv
        if(s == 1 & l == 1 & k == 1){
          write_csv(c_df,file=paste("Result/p03/Xpc-all/Xpc-pcd-all-",d,"-",n,".csv",sep = ""),
                    append = TRUE,col_names = TRUE)  
        }else{
          write_csv(c_df,file=paste("Result/p03/Xpc-all/Xpc-pcd-all-",d,"-",n,".csv",sep = ""),
                    append = TRUE)  
        }
        
      }
    }
  }
  
  }
}



end=Sys.time()
end-start
