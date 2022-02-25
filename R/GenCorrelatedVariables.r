rm(list=ls())


nvectors <- 500000       # Number of datasets to simulate. 

CorrZ <- matrix(0,3,3)   # CorrZ will be the correlation matrix of the standard normal vector.
CorrZ[1,2] <- -0.29      
CorrZ[1,3] <- -0.21
CorrZ[2,3] <- -0.31
CorrZ <- CorrZ + t(CorrZ) 
diag(CorrZ)<-1

muZ <- c(0,0,0)  

library(mvtnorm)
Z <- rmvnorm(nvectors,muZ,CorrZ)

U <- pnorm(Z)

Y <- matrix(cbind(qnorm(U[,1]),qbinom(U[,2],1,0.2),qpois(U[,3],3)),nvectors,3)

colMeans(Y)

cor(Y[,1],Y[,2])
cor(Y[,1],Y[,3])
cor(Y[,2],Y[,3])

