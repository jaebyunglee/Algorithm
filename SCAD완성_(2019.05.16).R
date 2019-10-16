rm(list=ls())
set.seed(1)
n = 1000; p = 30
x.mat = matrix(rnorm(n*p),ncol=p); b.vec = 1/(1:p)
exb.vec = exp(drop(x.mat%*%b.vec)); p.vec = exb.vec/(1+exb.vec)
y.vec = rbinom(n,1,prob=p.vec)
iter.max = 1e+3;eps = 1e-05


###concave function
# lam = 1; b = -8; a = 4
# s.fun = function(lam,a,b){
# fa = lam*abs(b)
# fb = (a*lam*(abs(b)-lam)-(abs(b)^2-lam^2)/2)/(a-1)+lam^2
# fc = (a-1)*lam^2/2+lam^2
# ret = fa*I(0<=abs(b)&abs(b)<lam)+fb*I(lam<=abs(b)&abs(b)<=a*lam)+fc*I(abs(b)>a*lam)-lam*abs(b)
#   return(ret)
# }
# 
# 
# bb.vec = NULL
# l.vec = seq(-5,5,length.out = 1000)
# for(b in l.vec){
#   bb.vec = c(bb.vec,s.fun(lam,a,b))
# }
# plot(bb.vec)



###concave grad function
ccav.fun = function(b.vec,lam,a){
  ret = rep(0,length(b.vec))
  ret = sign(b.vec)*ifelse(a*lam-abs(b.vec)>0,a*lam-abs(b.vec),0)/(a-1)-sign(b.vec)*lam
  if(sum(abs(b.vec)<=lam)!=0){ret[abs(b.vec)<=lam] = 0}
  return(ret)
}

###loss
loss.fun = function(y.vec,x.mat,b.vec){
  xb.vec = drop(x.mat%*%b.vec); ret = -sum(y.vec*xb.vec)+sum(log(1+exp(xb.vec)))  
  return(ret)
}
###scad.loss
scad.loss = function(y.vec,x.mat,b.vec,lam,a){
  fa = lam*abs(b.vec)
  fb = (a*lam*(abs(b.vec)-lam)-(abs(b.vec)^2-lam^2)/2)/(a-1)+lam^2
  fc = (a-1)*lam^2/2+lam^2
  ret = loss.fun(y.vec,x.mat,b.vec) + sum(fa*I(0<=abs(b.vec)&abs(b.vec)<lam)+fb*I(lam<=abs(b.vec)&abs(b.vec)<=a*lam)+fc*I(abs(b.vec)>a*lam)-lam*abs(b.vec))
  return(ret)
}

###grad
grad.fun = function(y.vec,x.mat,b.vec){
  xb.vec=drop(x.mat%*%b.vec);xb.vec = pmin(xb.vec,100)
  exb.vec = exp(xb.vec); p.vec = exb.vec/(1+exb.vec)
  ret = -drop(t(x.mat)%*%(y.vec-p.vec))/length(y.vec)
  return(ret)
}
###hess
hess.fun = function(y.vec,x.mat,b.vec){
  xb.vec = drop(x.mat%*%b.vec);xb.vec = pmin(xb.vec,100)
  exb.vec = exp(xb.vec); p.vec = exb.vec/(1+exb.vec)
  ret = t(x.mat)%*%diag(p.vec*(1-p.vec))%*%x.mat/length(y.vec) 
  return(ret)
}

###kkt condition
kkt.fun = function(g.vec,b.vec,lam,eps){
  kkt1 = prod(abs(g.vec[b.vec!=0]+lam*sign(b.vec[b.vec!=0]))<eps)
  kkt2 = prod(abs(g.vec[b.vec==0])-lam<eps)
  return(kkt1&kkt2)
}


###quadratic approximation
quad.lasso.fun = function(q.mat,l.vec,b.vec,lam,iter.max=1e+02,eps=1e-07){
  f.vec = rep(NA,iter.max)
  for(iter in 1:iter.max){
    for(j in 1:length(b.vec)){
      a = 2*q.mat[j,j]; b = -2*sum(q.mat[j,-j]*b.vec[-j])-l.vec[j]
      #print(cbind(lam,b))
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

###lasso
lasso.fun = function(y.vec,x.mat,b.vec,lam,iter.max=1e+3,eps=1e-5){
  for(iter in 1:iter.max){#iter
    h.mat = hess.fun(y.vec,x.mat,b.vec)
    g.vec = grad.fun(y.vec,x.mat,b.vec)
    q.mat = h.mat/2
    l.vec = g.vec-drop(h.mat%*%b.vec)
    b.vec = quad.lasso.fun(q.mat,l.vec,b.vec,lam,iter.max,eps)$b.vec
    g.vec = grad.fun(y.vec,x.mat,b.vec)
    conv = kkt.fun(g.vec,b.vec,lam,eps)
    if(conv) break
  }#iter
  #print(conv)
  #print(g.vec)
  return(list(b.vec=b.vec,g.vec=g.vec,lam=lam))
}

### active set control - lasso
act.lasso.fun = function(y.vec,x.mat,lam.vec,iter.max=1e+3,eps=1e-5){
  b.mat = NULL ;  b.vec = rep(0,ncol(x.mat))
  g.vec = grad.fun(y.vec,x.mat,b.vec)
  for(lam in lam.vec){#lam
    for(iter in 1:iter.max){#iter
      set1 = b.vec!=0
      set2 = abs(g.vec)>=(lam-eps)
      set = set1|set2
      if(sum(set)==0){set = 1:length(b.vec)}
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

### SCAD
scad.fun = function(y.vec,x.mat,b.vec,lam,a,iter.max=1e+3,eps=1e-5){
  for(c.id in 1:iter.max){#cccp
    #print(c.id)
    fb.val = scad.loss(y.vec,x.mat,b.vec,lam,a)
    #print(scad.loss(y.vec,x.mat,b.vec,lam,a))
    g.vec = grad.fun(y.vec,x.mat,b.vec)
    c.vec = ccav.fun(b.vec,lam,a)
    conv = kkt.fun(g.vec+c.vec,b.vec,lam,eps)
    if(conv) break 
    u.vec = b.vec*0
    for(q.id in 1:iter.max){#quad
      #print(q.id)
      h.mat = hess.fun(y.vec,x.mat,u.vec) 
      g.vec = grad.fun(y.vec,x.mat,u.vec)
      q.mat = h.mat/2
      l.vec = g.vec-drop(h.mat%*%u.vec)+c.vec
      u.vec = quad.lasso.fun(q.mat,l.vec,u.vec,lam,iter.max,eps)$b.vec
      g.vec = grad.fun(y.vec,x.mat,u.vec)
      conv = kkt.fun(g.vec+c.vec,u.vec,lam,eps)
      if(conv) break 
    }#quad
    fu.val = scad.loss(y.vec,x.mat,u.vec,lam,a)
    if(fu.val>fb.val){
      for(h in seq(0,0.01,length.out = 100)){
        b.vec = h*u.vec+(1-h)*b.vec
        # print(scad.loss(y.vec,x.mat,u.vec,lam,a))
        # print(b.vec[b.vec!=0])
      }
    } else { b.vec = u.vec }
  }#cccp
  return(list(b.vec=b.vec,"g.vec+c.vec"=g.vec+c.vec,lam=lam,a=a,conv=conv))
}


### active set control - SCAD
act.scad.fun = function(y.vec,x.mat,lam.vec,a,iter.max=1e+3,eps=1e-5){
  b.mat = NULL ;  b.vec = rep(0,ncol(x.mat))
  g.vec = grad.fun(y.vec,x.mat,b.vec)
  for(lam in lam.vec){#lam
    #print(lam)
    c.vec = ccav.fun(b.vec,lam,a)
    for(iter in 1:iter.max){#iter
      set1 = b.vec!=0
      set2 = abs(g.vec+c.vec)>=(lam-eps)
      set = set1|set2
      if(sum(set)==0){set = 1:length(b.vec)}
      ax.mat = x.mat[,set,drop=F]
      ab.vec = b.vec[set]
      fit = scad.fun(y.vec,ax.mat,ab.vec,lam,a)
      b.vec[set] = fit$b.vec
      g.vec = grad.fun(y.vec,x.mat,b.vec)
      c.vec = ccav.fun(b.vec,lam,a)
      conv = kkt.fun(g.vec+c.vec,b.vec,lam,eps)
      if(conv) break
    }#iter
    b.mat = cbind(b.mat,b.vec)
  }#lam
  return(list(b.mat=b.mat,lam.vec=lam.vec))
}
lam.vec = 10:1
act.lasso.fun(y.vec,x.mat,lam.vec=lam.vec,iter.max,eps)
### check
 lam.vec = seq(0.1,0.01,length.out=11)
# scad.fun(y.vec,x.mat,rep(0,p),lam=lam.vec[11],a=3,iter.max=1e+3,eps=1e-5)$b.vec
# act.scad.fun(y.vec,x.mat,lam.vec,a=3,iter.max=1e+3,eps=1e-5)$b.mat[,11]
# coef(glm(y.vec~.-1,data=data.frame(cbind(y.vec,x.mat)),family="binomial"))


scad1.fun = function(y.vec,x.mat,b.vec,lam,a,iter.max=1e+3,eps=1e-5){
  g.vec = grad.fun(y.vec,x.mat,b.vec)
  c.vec = ccav.fun(b.vec,lam,a)
  for(iter in 1:iter.max){#iter
    #print(iter)
    #print(scad.loss(y.vec,x.mat,b.vec,lam,a))
    h.mat = hess.fun(y.vec,x.mat,b.vec)
    q.mat = h.mat/2
    l.vec = g.vec-drop(h.mat%*%b.vec) + c.vec
    b.vec = quad.lasso.fun(q.mat,l.vec,b.vec,lam,iter.max,eps)$b.vec
    g.vec = grad.fun(y.vec,x.mat,b.vec)
    c.vec = ccav.fun(b.vec,lam,a)
    conv = kkt.fun(g.vec+c.vec,b.vec,lam,eps)
    if(conv) break
  }#iter
  return(list(b.vec=b.vec,g.vec=g.vec+c.vec,lam=lam))
}




scad.fun(y.vec,x.mat,rep(0,p),lam=0.01,a=3,iter.max=1e+3,eps=1e-5)
scad1.fun(y.vec,x.mat,rep(0,p),lam=0.01,a=3,iter.max=1e+3,eps=1e-5)
act.scad.fun(y.vec,x.mat,lam.vec,3,iter.max=1e+3,eps=1e-5)
coef(glm(y.vec~.-1,data = data.frame(cbind(y.vec,x.mat)),family="binomial"))




