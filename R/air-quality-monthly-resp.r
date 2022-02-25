rm(list=ls())

# input data
regTaiwan_M <- read.csv("DataExample/taiwan-monthly.csv",header = T)
regTaiwan_M <- data.frame(regTaiwan_M)

# ARRANGE THE DATA AS A LIST OF DATA SETS
regionsTaiwan_M <- as.character(unique(regTaiwan_M$city))
dlistTaiwan_M <- lapply(regionsTaiwan_M,function(x) regTaiwan_M[regTaiwan_M$city==x,])
names(dlistTaiwan_M) <- regionsTaiwan_M

## Kaohsing 

## Phase I  36 samples (monthly) of 30( or 31) observations
## Consider X1 = Ozone ; X2 = PM10_quality ; X3 = all death

##----------------- x-chart --------------------------------------------
data_M <- dlistTaiwan_M[[1]]

n <- 36

## Phase I 
x1 <- data_M$Ozone[1:n]

mu0 <- mean(x1)
sig0 <- 19.91

pval.x0 <- numeric()

for(i in 1:n){
  lower.x <- pnorm((x1[i]-mu0)/(sig0/sqrt(n)))
  upper.x <- ( 1- pnorm((x1[i]-mu0)/(sig0/sqrt(n))) )
  
  pval.x0[i] <- round(2*min(c(lower.x,upper.x)),5)
}

## Phase II 

z <- data_M$Ozone[(n+1):(2*n)] 

pval.x1 <- numeric()
for(i in 1:n){
  lower.x <- pnorm((z[i]-mu0)/(sig0/sqrt(n)))
  upper.x <- ( 1- pnorm((z[i]-mu0)/(sig0/sqrt(n))) )
  
  pval.x1[i] <- round(2*min(c(lower.x,upper.x)),5)
}

##----------------- p-chart --------------------------------------------

## Phase I
zp0 <- data_M$PM10[1:n]

p0 <- mean(zp0)


pval.p0 <- numeric()
for(i in 1:n){
  lower.p <- pnorm((zp0[i]-p0)/sqrt(p0*(1-p0)/n))
  upper.p <- ( 1- pnorm((zp0[i]-p0)/sqrt(p0*(1-p0)/n)) )
  
  pval.p0[i] <- round(2*min(c(lower.p,upper.p)),5)
}

## Phase II
zp1 <- data_M$PM10[(n+1):(2*n)]

pval.p1 <- numeric()
for(i in 1:n){
  lower.p <- pnorm((zp1[i]-p0)/sqrt(p0*(1-p0)/n))
  upper.p <- ( 1- pnorm((zp1[i]-p0)/sqrt(p0*(1-p0)/n)) )
  
  pval.p1[i] <- round(2*min(c(lower.p,upper.p)),5)
}


##----------------- c-chart --------------------------------------------

## Phase I
zc0 <- data_M$resp[1:n]
lam0 <- mean(zc0)

pval.c0 <- numeric()
for(i in 1:n){
  lower.c <- pnorm((zc0[i]-lam0)/sqrt(lam0/n))
  upper.c <- ( 1- pnorm((zc0[i]-lam0)/sqrt(lam0/n)) )
  
  pval.c0[i] <- round(2*min(c(lower.c,upper.c)),5) 
}


## Phase II 

zc1 <- data_M$resp[(n+1):(2*n)]


pval.c1 <- numeric()
for(i in 1:n){
  lower.c <- pnorm((zc1[i]-lam0)/sqrt(lam0/n))
  upper.c <- ( 1- pnorm((zc1[i]-lam0)/sqrt(lam0/n)) )
  
  pval.c1[i] <- round(2*min(c(lower.c,upper.c)),5) 
}


##------------- Combined ------------------------------------


##------- Plot ---------------------- 

A <- cbind(pval.x0,pval.p0,pval.c0)
A

A_ord <- matrix(0,n,3)
for(i in 1:n){
  A_ord[i,] <- sort(A[i,])
}


A1 <- cbind(pval.x1,pval.p1,pval.c1)
A1

A1_ord <- matrix(0,n,3)
for(i in 1:n){
  A1_ord[i,] <- sort(A1[i,])
}

B1 <- rbind(A,A1)
B <- rbind(A_ord,A1_ord)


## Combine three chart on the same figure
par(mai=c(1.0,1.2,0.3,0.3))
x <- c(1,6,12,18,24,30,36,42,48,54,60,66,72)
y <- c(0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1)
y1 <- log(y*500000)
plot(log(B1[,1]*500000),type = "b",pch=16, xaxt = "n", yaxt = "n",
     xlim = c(2.5,71),ylim = c(1,13.5),col="red",lty=1,lwd=1,
     xlab = "Sample Sequence",ylab = expression(paste("P-value")),cex.lab=1.8,cex.axis=1.3)
points(log(B1[,2]*500000),col="blue",type = "b",pch=16,lty=1,lwd=1.3)
points(log(B1[,3]*500000),type = "b",pch=16,lty=1,lwd=1.6)
legend(3.2,5,legend = c(expression(X[1]),expression(X[2]),expression(X[3])),
       col=c("red","blue",1),lty = c(1,1,1),lwd = c(1,1.3,1.6),cex = 1.3 )
axis(1, labels = x, at = x,cex.axis=1.1)
axis(2, labels = y, at = y1)
abline(h=log(0.005*500000),col=6)
abline(h=log(0.005/2*500000),col="chartreuse3")
abline(h=log(0.005/3*500000),col="steelblue")
abline(v=36)
text(31,1.15,"Phase I",cex = 3)
text(41,1.15,"Phase II",cex = 3)
text(3.2,8,expression(paste(LCL[(3)]," = 0.005")),cex = 1)
text(3.2,7.3,expression(paste(LCL[(2)]," = 0.0025")),cex = 1)
text(3.2,6.9,expression(paste(LCL[(1)]," = 0.0016")),cex = 1)

rej_c <- which(B1[,3] < 0.00166)
points(rej_c,log(B1[rej_c,3]*500000),pch=8,cex=1.2)
text(rej_c,log(B1[rej_c,3]*500000)-0.2,rej_c,cex=1.1)



## Bonferroni's adjustment 
par(mai=c(1.0,1.2,0.3,0.3))
x <- c(1,6,12,18,24,30,36,42,48,54,60,66,72)
y <- c(0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1)
y1 <- log(y*500000)
plot(log(B1[,1]*500000),type = "b",pch=16, xaxt = "n", yaxt = "n",
     xlim = c(2.5,71),ylim = c(1,13.5),col="red",lty=1,lwd=1,
     xlab = "Sample Sequence",ylab = expression(paste("P-value")),cex.lab=1.8,cex.axis=1.3)
points(log(B1[,2]*500000),col="blue",type = "b",pch=16,lty=1,lwd=1.3)
points(log(B1[,3]*500000),type = "b",pch=16,lty=1,lwd=1.6)
legend(3.2,5,legend = c(expression(X[1]),expression(X[2]),expression(X[3])),
       col=c("red","blue",1),lty = c(1,1,1),lwd = c(1,1.3,1.6),cex = 1.3 )
axis(1, labels = x, at = x,cex.axis=1.1)
axis(2, labels = y, at = y1)
abline(h=log(0.005/3*500000),col="steelblue")
abline(v=36)
text(31,1.15,"Phase I",cex = 3)
text(41,1.15,"Phase II",cex = 3)
text(3.2,7,"LCL = 0.0016",cex = 1.2)
rej_c <- which(B1[,3] < 0.00166)
points(rej_c,log(B1[rej_c,3]*500000),pch=8,cex=1.2)
text(rej_c,log(B1[rej_c,3]*500000)-0.2,rej_c,cex=1.1)






