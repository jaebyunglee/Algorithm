rm(list=ls())
library(MASS)
library(mvtnorm)
set.seed(1)
n = 5000
x1.mat = mvrnorm(n, mu = c(1,2.25),Sigma = matrix(c(1.75,1,1,1.75),2,2) )
x2.mat = mvrnorm(n, mu = c(12.9,10.3),Sigma = matrix(c(3.2,1.25,1.25,3.2),2,2))
I=rbinom(n,size=1,prob = 0.3)
x.mat = I*x1.mat + (1-I)*x2.mat
plot(x.mat)


t = 0.4
u1 = c(0,0)
u2 = c(1,1)
Sigma1 = matrix(c(1,0.8,0.8,1),2,2)
Sigma2 = matrix(c(1.5,0.3,0.3,1.5),2,2)
#############################################



for(i in 1:100){
  d1=dmvnorm(x.mat,mean = u1,sigma = Sigma1)
  d2=dmvnorm(x.mat,mean = u2,sigma = Sigma2)
  #E-step
  T1=(t*d1)/(t*d1+(1-t)*d2)
  T2=(1-t)*d2/(t*d1+(1-t)*d2)
  
  #M-step
  nu1 = colSums(T1*x.mat)/sum(T1)
  nu2 = colSums(T2*x.mat)/sum(T2)
  s1 = scale(x.mat,center = nu1,scale = FALSE)
  s2 = scale(x.mat,center = nu2,scale = FALSE)
  nSigma1=(t(s1)%*%diag(T1)%*%(s1))/sum(T1)
  nSigma2=(t(s2)%*%diag(T2)%*%(s2))/sum(T2)
  nt = sum(T1)/n

  if(sum(abs(u1-nu1))+sum(abs(u2-nu2))+sum(abs(Sigma1-nSigma1))+sum(abs(Sigma2-nSigma2))<1e-10) break
  print(c(i,sum(abs(u1-nu1))+sum(abs(u2-nu2))+sum(abs(Sigma1-nSigma1))+sum(abs(Sigma2-nSigma2))))
  cat('accracy',sum(diag(table(T1>0.5,I)))/sum(table(T1>0.5,I)),'\n')
  t = nt ; u1 = nu1 ; u2 = nu2 ; Sigma1 = nSigma1 ; Sigma2 = nSigma2
}
u1; u2
Sigma1; Sigma2


par(mfrow=c(1,2))
plot(x.mat)
plot(x.mat[T1>0.5,],ylim = c(-3,16),xlim=c(-3,20),col="red")
points(x.mat[T1<=0.5,],col="blue")
