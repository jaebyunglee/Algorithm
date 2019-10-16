rm(list=ls())
set.seed(1)
n = 100; p = 10
x.mat = matrix(rnorm(n*p),ncol=p); b.vec = 1/(1:p)
exb.vec = exp(drop(x.mat%*%b.vec)); p.vec = exb.vec/(1+exb.vec)
y.vec = rbinom(n,1,prob=p.vec)

###loss
loss.fun = function(y.vec,x.mat,b.vec){
  xb.vec = drop(x.mat%*%b.vec); ret = -sum(y.vec*xb.vec)+sum(log(1+exp(xb.vec)))  
  return(ret)
}
###gradient
grad.fun = function(y.vec,x.mat,b.vec){
  exb.vec = exp(drop(x.mat%*%b.vec)); p.vec = exb.vec/(1+exb.vec)
  ret = -drop(t(x.mat)%*%(y.vec-p.vec))
  return(ret)
}
###hessian
hess.fun = function(y.vec,x.mat,b.vec){
  exb.vec = exp(drop(x.mat%*%b.vec)); p.vec = exb.vec/(1+exb.vec)
  ret = t(x.mat)%*%diag(p.vec*(1-p.vec))%*%x.mat 
  return(ret)
}

eps = 1e-10; iter.max = 100
b.vec = rep(0,p); lam = 1
q.mat = hess.fun(y.vec,x.mat,b.vec)
l.vec = grad.fun(y.vec,x.mat,b.vec)

###quad.lasso.fun & kkt condition
quad.lasso.fun = function(q.mat,l.vec,lam,b.vec,iter.max,eps){
  f.vec = rep(NA,iter.max)
  for(iter in 1:iter.max){
    for(j in 1:length(b.vec)){
      a = 2*q.mat[j,j]; b = -2*sum(q.mat[j,-j]*b.vec[-j])-l.vec[j]
      if(abs(b)<lam){ b.vec[j] = 0 
      } else { b.vec[j] = sign(b)*(abs(b)-lam)/a } 
    }
    f.vec[iter] = t(b.vec)%*%q.mat%*%b.vec+sum(l.vec*b.vec)+lam*sum(abs(b.vec))
    grad = 2*q.mat%*%b.vec+l.vec
    kkt1 = prod((grad[b.vec!=0] + lam*sign(b.vec[b.vec!=0]))<eps)
    kkt2 = prod(abs(grad[b.vec==0])-lam<eps)
    print(c(kkt1,kkt2))
    if(kkt1&kkt2) break 
  }
  return(list(b.vec=b.vec,f.vec=f.vec[1:iter]))
}

quad.lasso.fun(q.mat,l.vec,lam,b.vec,iter.max,eps)


###home work
grad = 2*q.mat%*%b.vec+l.vec
qq.mat = matrix(q.mat[1,1],1,1)
ll.vec = l.vec[1]
fit = quad.lasso.fun(qq.mat,ll.vec,lam,b.vec[1],iter.max,eps)
b.vec[1] = fit$b.vec 
grad = 2*q.mat%*%b.vec+l.vec

qq.mat = matrix(q.mat[1:2,1:2],2,2)
ll.vec = l.vec[1:2]
fit = quad.lasso.fun(qq.mat,ll.vec,lam,b.vec[1:2],iter.max,eps)
b.vec[1:2] = fit$b.vec 
grad = 2*q.mat%*%b.vec+l.vec
