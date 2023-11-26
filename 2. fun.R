

spike_n_slab <- function(y, X, n.post = 500, n.disp = 50, pi0 = 0.01, tau2 = 100, init = NULL){
  
  Z.mat = matrix(0, n.post, p)
  beta.mat = matrix(0, n.post, p)
  sig2.mat = matrix(0, n.post, 1)
  
  # hyperparameters
  # pi0 = 0.01
  # pi0 = 0.1
  # pi0 = 0.25
  # pi0 = 0.5
  # tau2 = 100
  
  # initial values
  if(is.null(init)){
    Z.vec = c(1,1,1, rep(0,p-3))
    beta.vec = c(1,1,1, rep(0,p-3))
    sigma2.val = 2
  }
  
  for(iter in 1:n.post){
    
    # sample Z
    for(j in 1:p){
      Z.vec[j] = 1
      Z.vec.mj = Z.vec
      Z.vec.mj[j] = 0
      
      X_Z = X[, which(Z.vec==1), drop=FALSE]
      gram.X_Z = t(X_Z) %*% X_Z
      det.Z.term = determinant(tau2*gram.X_Z + diag(sum(Z.vec)), logarithm = TRUE)$mod
      hat.Z.term = X_Z%*%solve(gram.X_Z + 1/tau2*diag(sum(Z.vec)) )%*%t(X_Z)
      
      if(sum(Z.vec.mj)==0){
        det.Z.mj.term = 0
        hat.Z.mj.term = -diag(n)
      }else{
        X_Z_mj = X[, which(Z.vec.mj==1), drop=FALSE]
        gram.X_Z_mj = t(X_Z_mj) %*% X_Z_mj
        det.Z.mj.term = determinant(tau2*gram.X_Z_mj + diag(sum(Z.vec.mj)), logarithm = TRUE)$mod
        hat.Z.mj.term = X_Z_mj%*%solve(gram.X_Z_mj + 1/tau2*diag(sum(Z.vec.mj)) )%*%t(X_Z_mj)
      }
      
      log.qj = log(pi0) - log(1-pi0) + 0.5*(det.Z.mj.term - det.Z.term) +
        0.5*1/sigma2.val*t(y)%*%(hat.Z.term - hat.Z.mj.term)%*%y
      qj = exp(log.qj[1])
      if(qj==Inf){
        Z.vec[j] = 1
      }else{
        Z.vec[j] = rbinom(1, size=1, prob=qj/(1+qj))
      }
    }
    
    # sample beta
    if(sum(Z.vec)==0){
      beta.vec = rep(0, p)
    }else{
      X_Z = X[, which(Z.vec==1), drop=FALSE]
      gram.X_Z = t(X_Z) %*% X_Z
      inv.Sg.k = gram.X_Z + 1/tau2*diag(sum(Z.vec))
      Sg.k = solve(inv.Sg.k)
      Mu.k = Sg.k %*% t(X_Z)%*%y
      beta.vec[which(Z.vec==1)] = rmvnorm(1, mean=Mu.k, sigma=sigma2.val*Sg.k)
      beta.vec[which(Z.vec==0)] = 0
    }
    
    # sample sigma2
    if(sum(Z.vec)==0){
      sigma2.val = rinvgamma(n=1, 
                             shape=n/2, 
                             rate=sum((y)^2)/2)
    }else{
      sigma2.val = rinvgamma(n=1, 
                             shape=(n+sum(Z.vec))/2, 
                             rate=( sum((y- X_Z%*%beta.vec[which(Z.vec==1)])^2) + 1/tau2*sum((beta.vec)^2) )/2)
    }
    
    Z.mat[iter,] = Z.vec
    beta.mat[iter,] = beta.vec
    sig2.mat[iter,] = sigma2.val
    
    if(iter %% n.disp == 0) print(iter)
  }
  
  res = list(Z.mat=Z.mat, beta.mat=beta.mat, sig2.mat=sig2.mat)
  return(res)
}






