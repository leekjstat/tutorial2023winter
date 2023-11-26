
library(mvtnorm)
library(invgamma)
library(ggplot2)
library(ggpubr) # for ggarrange
source("fun.R")

n = 50
p = 100

sigma2.0 = 1

set.seed(1)
X = matrix(rnorm(n*p), n, p)
eps = matrix(rnorm(n, sd=sqrt(sigma2.0)), n, 1)

Z.0 = c(1,1,1, rep(0,p-3))
beta.0 = rep(0, p)
beta.0[which(Z.0==1)] = 1

y = X%*%beta.0 + eps
var(X%*%beta.0)/sigma2.0 # SNR


# PIP for different pi0's
res1 = spike_n_slab(y, X, n.post = 500, n.disp = 100, pi0 = 0.01, tau2 = 100, init = NULL)
res2 = spike_n_slab(y, X, n.post = 500, n.disp = 100, pi0 = 0.2, tau2 = 100, init = NULL)
res3 = spike_n_slab(y, X, n.post = 500, n.disp = 100, pi0 = 0.4, tau2 = 100, init = NULL)
res4 = spike_n_slab(y, X, n.post = 500, n.disp = 100, pi0 = 0.6, tau2 = 100, init = NULL)

PIP1 = colMeans(res1$Z.mat)
PIP2 = colMeans(res2$Z.mat)
PIP3 = colMeans(res3$Z.mat)
PIP4 = colMeans(res4$Z.mat)

df1 = data.frame(x=1:p, PIP=PIP1, ind=as.factor(1-Z.0))
df2 = data.frame(x=1:p, PIP=PIP2, ind=as.factor(1-Z.0))
df3 = data.frame(x=1:p, PIP=PIP3, ind=as.factor(1-Z.0))
df4 = data.frame(x=1:p, PIP=PIP4, ind=as.factor(1-Z.0))


gp1 = ggplot(df1,aes(x=x,y=PIP,col=ind)) + geom_point() + theme(legend.position='none') + 
  labs(title=expression(pi[0]==0.01), y="Posterior inclusion prob.") + theme(plot.title = element_text(hjust = 0.5, size=15))
gp2 = ggplot(df2,aes(x=x,y=PIP,col=ind)) + geom_point() + theme(legend.position='none') + 
  labs(title=expression(pi[0]==0.2), y="Posterior inclusion prob.") + theme(plot.title = element_text(hjust = 0.5, size=15))
gp3 = ggplot(df3,aes(x=x,y=PIP,col=ind)) + geom_point() + theme(legend.position='none') + 
  labs(title=expression(pi[0]==0.4), y="Posterior inclusion prob.") + theme(plot.title = element_text(hjust = 0.5, size=15))
gp4 = ggplot(df4,aes(x=x,y=PIP,col=ind)) + geom_point() + theme(legend.position='none') + 
  labs(title=expression(pi[0]==0.6), y="Posterior inclusion prob.") + theme(plot.title = element_text(hjust = 0.5, size=15))
ggarrange(gp1, gp2, gp3, gp4, ncol=2, nrow=2)



# PIP for different tau2's
res5 = spike_n_slab(y, X, n.post = 500, n.disp = 100, pi0 = 0.01, tau2 = 0.01, init = NULL)
res6 = spike_n_slab(y, X, n.post = 500, n.disp = 100, pi0 = 0.01, tau2 = 1, init = NULL)
res7 = spike_n_slab(y, X, n.post = 500, n.disp = 100, pi0 = 0.01, tau2 = 100, init = NULL)
res8 = spike_n_slab(y, X, n.post = 500, n.disp = 100, pi0 = 0.01, tau2 = 1000000, init = NULL)

PIP5 = colMeans(res5$Z.mat)
PIP6 = colMeans(res6$Z.mat)
PIP7 = colMeans(res7$Z.mat)
PIP8 = colMeans(res8$Z.mat)

df5 = data.frame(x=1:p, PIP=PIP5, ind=as.factor(1-Z.0))
df6 = data.frame(x=1:p, PIP=PIP6, ind=as.factor(1-Z.0))
df7 = data.frame(x=1:p, PIP=PIP7, ind=as.factor(1-Z.0))
df8 = data.frame(x=1:p, PIP=PIP8, ind=as.factor(1-Z.0))

gp5 = ggplot(df5,aes(x=x,y=PIP,col=ind)) + geom_point() + theme(legend.position='none') + 
  labs(title=expression(tau^2==0.01), y="Posterior inclusion prob.") + theme(plot.title = element_text(hjust = 0.5, size=15))
gp6 = ggplot(df6,aes(x=x,y=PIP,col=ind)) + geom_point() + theme(legend.position='none') + 
  labs(title=expression(tau^2==1), y="Posterior inclusion prob.") + theme(plot.title = element_text(hjust = 0.5, size=15))
gp7 = ggplot(df7,aes(x=x,y=PIP,col=ind)) + geom_point() + theme(legend.position='none') + 
  labs(title=expression(tau^2==100), y="Posterior inclusion prob.") + theme(plot.title = element_text(hjust = 0.5, size=15))
gp8 = ggplot(df8,aes(x=x,y=PIP,col=ind)) + geom_point() + theme(legend.position='none') + 
  labs(title=expression(tau^2==1000000), y="Posterior inclusion prob.") + theme(plot.title = element_text(hjust = 0.5, size=15))
ggarrange(gp5, gp6, gp7, gp8, ncol=2, nrow=2)


bm5 = colMeans(res5$beta.mat)
bm6 = colMeans(res6$beta.mat)
bm7 = colMeans(res7$beta.mat)
bm8 = colMeans(res8$beta.mat)
bm.df5 = data.frame(x=1:p, bm=bm5, ind=as.factor(1-Z.0))
bm.df6 = data.frame(x=1:p, bm=bm6, ind=as.factor(1-Z.0))
bm.df7 = data.frame(x=1:p, bm=bm7, ind=as.factor(1-Z.0))
bm.df8 = data.frame(x=1:p, bm=bm8, ind=as.factor(1-Z.0))
bm.gp5 = ggplot(bm.df5,aes(x=x,y=bm,col=ind)) + geom_point() + theme(legend.position='none') + 
  labs(title=expression(tau^2==0.01), y="Coefficient") + theme(plot.title = element_text(hjust = 0.5, size=15))
bm.gp6 = ggplot(bm.df6,aes(x=x,y=bm,col=ind)) + geom_point() + theme(legend.position='none') + 
  labs(title=expression(tau^2==1), y="Coefficient") + theme(plot.title = element_text(hjust = 0.5, size=15))
bm.gp7 = ggplot(bm.df7,aes(x=x,y=bm,col=ind)) + geom_point() + theme(legend.position='none') + 
  labs(title=expression(tau^2==100), y="Coefficient") + theme(plot.title = element_text(hjust = 0.5, size=15))
bm.gp8 = ggplot(bm.df8,aes(x=x,y=bm,col=ind)) + geom_point() + theme(legend.position='none') + 
  labs(title=expression(tau^2==1000000), y="Coefficient") + theme(plot.title = element_text(hjust = 0.5, size=15))
ggarrange(bm.gp5, bm.gp6, bm.gp7, bm.gp8, ncol=2, nrow=2)


est.Z5 = as.numeric(PIP5 > 0.5)
est.Z6 = as.numeric(PIP6 > 0.5)
est.Z7 = as.numeric(PIP7 > 0.5)
est.Z8 = as.numeric(PIP8 > 0.5)

if(sum(est.Z5)==0){
  
}else{
  mu5 = solve(t(X[,which(est.Z5==1)])%*%X[,which(est.Z5==1)]+1/0.01*diag(sum(est.Z5)) ) %*% t(X[,which(est.Z5==1)]) %*% y
}

mu6 = solve(t(X[,which(est.Z6==1)])%*%X[,which(est.Z6==1)]+1/1*diag(sum(est.Z6)) ) %*% t(X[,which(est.Z6==1)]) %*% y
mu7 = solve(t(X[,which(est.Z7==1)])%*%X[,which(est.Z7==1)]+1/100*diag(sum(est.Z7)) ) %*% t(X[,which(est.Z7==1)]) %*% y
mu8 = solve(t(X[,which(est.Z8==1)])%*%X[,which(est.Z8==1)]+1/1000000*diag(sum(est.Z8)) ) %*% t(X[,which(est.Z8==1)]) %*% y


