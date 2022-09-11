library(MASS)

# MH algorithm
MH <- function(f=fstar,s_matrix= sigma_matrix,x0=init,n=5e3,nburn=5e2) {
  Xres <- matrix(0, nrow = (n+nburn), ncol = 2)
  Xres[1,] <- x0
  accept <- 0
  for (i in 2:(n+nburn)){
    Y <- mvrnorm(n = 1, Xres[i-1,], s_matrix)
    if (runif(1) <= min(1, exp(log(f(Y)) -log(f(Xres[i-1,])) ))) {
      Xres[i,] <- Y
      accept <- accept+1
    }
    else
      Xres[i,] <- Xres[i-1,]
  }
  
  #percentage of accepted values
  print (accept/(n+nburn))
  return (Xres[-(1:nburn),])
}


#Adaptive Proposal algorithm

MH_AP <- function(f=fstar,s_matrix=sigma_matrix,x0=init,
                  n=1e3,h=200,U=200,c=2.4/sqrt(2)) {
  Xres <- matrix(0, nrow = (n), ncol = 2)
  Xres[1,] <- x0
  accept <- 0
  for (i in 2:U){
    #standard MH
    Y <- mvrnorm(n = 1, Xres[i-1,], s_matrix)
    if (runif(1) <= min(1, exp(log(f(Y)) -log(f(Xres[i-1,])) ))) {
      Xres[i,] <- Y
      accept <- accept+1
    }
    else
      Xres[i,] <- Xres[i-1,]
  }
  
  #K from the original paper
  K = Xres[(U-h+1):U,]
  for (j in (U+1):n){
    if (j%%U==1){
      K_tilde<-matrix(0,nrow=h,ncol=2)
      K = Xres[((j-1)-h+1):(j-1),]
      K_tilde[,1]=K[,1]-mean(K[,1])
      K_tilde[,2]=K[,2]-mean(K[,2])
      R<-1/(h-1)*t(K_tilde)%*%K_tilde
    }
    i = j
    Y <- mvrnorm(n = 1, Xres[i-1,], c^2*R)
    if (runif(1) <= min(1, exp(log(f(Y)) -log(f(Xres[i-1,]))))) {
      Xres[i,] <- Y
      accept <- accept+1
    }
    else
      Xres[i,] <- Xres[i-1,]
  }
  #percentage of accepted values
  print (accept/(n))
  return (Xres)
}

fstar=function(X) exp(-X[1]^2/100 - (X[2]+ 3*X[1]^2/100 - 3)^2 )

sigma_matrix=matrix(c(50,1,1,2.5),nrow=2,ncol=2)
init=c(0,0)


Samples1<-MH(n=1e5, x0=c(2.5,2.5),s_matrix=matrix(c(7.348^2,-1.2,-1.2,2.097^2),nrow=2,ncol=2))
Samples<-MH_AP(n = 1e5)


mc_chains_x_classic = mcmc(Samples1[,1])
ESS_x_classic = effectiveSize(mc_chains_x_classic)

mc_chains_y_classic = mcmc(Samples1[,2])
ESS_y_classic = effectiveSize(mc_chains_y_classic)
ESS

mc_chains_x_adaptive = mcmc(Samples[,1])
ESS_x_adaptive = effectiveSize(mc_chains_x_adaptive)

mc_chains_y_adaptive = mcmc(Samples1[,2])
ESS_y_adaptive = effectiveSize(mc_chains_y_adaptive)



