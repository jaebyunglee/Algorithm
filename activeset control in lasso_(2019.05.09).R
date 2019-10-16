rm(list=ls())

set.seed(1)
n = 100; p = 10
x.mat = matrix(rnorm(n*p),ncol=p); b.vec = 1/(1:p)
exb.vec = exp(drop(x.mat%*%b.vec)); p.vec = exb.vec/(1+exb.vec)
y.vec = rbinom(n,1,prob=p.vec)
iter.max = 1e+3;eps = 1e-05

###loss
loss.fun = function(y.vec,x.mat,b.vec){
  xb.vec = drop(x.mat%*%b.vec); ret = -sum(y.vec*xb.vec)+sum(log(1+exp(xb.vec)))  
  return(ret)
}
###grad
grad.fun = function(y.vec,x.mat,b.vec){
  exb.vec = exp(drop(x.mat%*%b.vec)); p.vec = exb.vec/(1+exb.vec)
  ret = -drop(t(x.mat)%*%(y.vec-p.vec))
  return(ret)
}
###hess
hess.fun = function(y.vec,x.mat,b.vec){
  exb.vec = exp(drop(x.mat%*%b.vec)); p.vec = exb.vec/(1+exb.vec)
  ret = t(x.mat)%*%diag(p.vec*(1-p.vec))%*%x.mat 
  return(ret)
}

###kkt condition
kkt.fun = function(g.vec,b.vec,lam,eps){
  kkt1 = prod(abs(g.vec[b.vec!=0] + lam*sign(b.vec[b.vec!=0]))<eps)
  kkt2 = prod(abs(g.vec[b.vec==0])-lam<eps)
  return(kkt1&kkt2)
}

###quadratic approximation
quad.lasso.fun = function(q.mat,l.vec,lam,b.vec,iter.max=1e+02,eps=1e-07){
  f.vec = rep(NA,iter.max)
  for(iter in 1:iter.max){
    for(j in 1:length(b.vec)){
      a = 2*q.mat[j,j]; b = -2*sum(q.mat[j,-j]*b.vec[-j])-l.vec[j]
      if(abs(b)<lam){ b.vec[j] = 0 
      } else { b.vec[j] = sign(b)*(abs(b)-lam)/a } 
    }
    f.vec[iter] = t(b.vec)%*%q.mat%*%b.vec+sum(l.vec*b.vec)+lam*sum(abs(b.vec))
    g.vec = 2*q.mat%*%b.vec+l.vec
    conv = kkt.fun(g.vec,b.vec,lam,eps)
    if(conv) break 
  }
  return(list(b.vec=b.vec,g.vec=g.vec,conv=conv,lam=lam))
}


###quadratic approximation lasso  ###SCAD homework
                                 #initial beta
lasso.fun = function(y.vec,x.mat,b.vec,lam,iter.max=1e+3,eps=1e-5){
  for(iter in 1:iter.max){#iter
    h.mat = hess.fun(y.vec,x.mat,b.vec)
    g.vec = grad.fun(y.vec,x.mat,b.vec)
    q.mat = h.mat/2
    l.vec = g.vec-drop(h.mat%*%b.vec)
    b.vec = quad.lasso.fun(q.mat,l.vec,lam,b.vec,iter.max,eps)$b.vec
    g.vec = grad.fun(y.vec,x.mat,b.vec)
    conv = kkt.fun(g.vec,b.vec,lam,eps)
    if(conv) break
  }#iter
  #print(conv)
  #print(g.vec)
 return(list(b.vec=b.vec,g.vec=g.vec,lam=lam))
}


###active set control
lam.vec = 10:1
act.lasso.fun = function(y.vec,x.mat,lam.vec,iter.max=1e+3,eps=1e-5){
  b.mat = NULL ;  b.vec = rep(0,ncol(x.mat))
  g.vec = grad.fun(y.vec,x.mat,b.vec)
  for(lam in lam.vec){#lam
    for(iter in 1:iter.max){#iter
      set1 = b.vec!=0
      set2 = abs(g.vec)>=(lam-eps)
      set = set1|set2
      ax.mat = x.mat[,set,drop=F]
      ab.vec = b.vec[set]
      fit = lasso.fun(y.vec,ax.mat,ab.vec,lam)
      b.vec[set] = fit$b.vec
      g.vec = grad.fun(y.vec,x.mat,b.vec)
      conv = kkt.fun(g.vec,b.vec,lam,eps)
      if(conv) break
    }#iter
    b.mat = cbind(b.mat,b.vec)
  }#lam
  return(list(b.mat=b.mat,lam.vec=lam.vec))
}


act.lasso.fun(y.vec,x.mat,lam.vec=10:1,iter.max,eps)
lasso.fun(y.vec,x.mat,10,b.vec=rep(0,p),iter.max,eps)
