#Parallel Tempering
g <- function(x) 0.4*exp(-(3-x)^2)+0.6*exp(-(30-x)^2)

#tempered version of original function
g_tempered<-function(x,t){
  exp(log(0.4*exp(-(3-x)^2)+0.6*exp(-(30-x)^2))/t)
}
q_exchange = function(a,b,m){
  if (a==1 & b ==2){
    return (1)
  }
  else if (a==m & b==(m-1)){
    return (1)
  }
  else if (abs(a-b)==1){
    return(0.5)
  }
  else {return (0)}
}
MH_Parallel_Tempering <- function(ftemp = g_tempered,sigma = c(0.5,2.8,1,1.5,4,8),x0 = c(15,25,0,10,
                                                                       45,5), n=1e4, temperatures = c(1,2,5,10,20,50),c=1) {
  Temps <- temperatures
  m <- length(Temps) #number of chains
  swap = matrix(0,nrow=1,ncol=m)
  
  Samples = list()
  for (i in 1:m){
    Samples[[i]] = matrix(0, nrow = n, ncol = 1)
    Samples[[i]][1,1] = x0[i]
  }
  count <- 0
  for (i in 2:n){
    #MH step
    for(j in 1:m){
      X <- Samples[[j]]
      Y = X[i-1,1] + rnorm(1,0,sigma[j])
      temp<-Temps[j]
      numerator = log(ftemp(Y,temp))
      denominator = log(ftemp(X[i-1,],temp))
      if (runif(1) <= min(1, exp(numerator-denominator))) {
        X[i,1] <- Y
        count <- count+1
      }
      else {
        X[i,1] <- X[i-1,1]
      }
      Samples[[j]] <- X
    }
    #Start Exchanging
    for (j in 1:m){
      if (j==1){
        #Follow the rule below
        #l = j
        #m = j+1
        fraction_1 = log(g_tempered(Samples[[j+1]][i,1],1))-log((g_tempered(Samples[[j]][i,1],1)))
        alpha = min(1,exp(fraction_1*(-1/Temps[j+1]+1/Temps[j]))*
                      q_exchange(j+1,j,m)/q_exchange(j,j+1,m))
        if (runif(1)<alpha*c){
          Exchange_Value = Samples[[j]][i,1]
          Samples[[j]][i,1]=Samples[[j+1]][i,1]
          Samples[[j+1]][i,1]=Exchange_Value
          swap[1,j] = swap[1,j] + 1
        }
      }
      else if (j==m){
        #Follow the rule below
        #l = j
        #m = j-1
        fraction_1 = log(g_tempered(Samples[[j-1]][i,1],1))-log((g_tempered(Samples[[j]][i,1],1)))
        alpha = min(1,exp(fraction_1*(-1/Temps[j-1]+1/Temps[j]))*
                      q_exchange(j-1,j,m)/q_exchange(j,j-1,m))
        if (runif(1)<alpha*c){
          Exchange_Value = Samples[[j]][i,1]
          Samples[[j]][i,1]=Samples[[j-1]][i,1]
          Samples[[j-1]][i,1]=Exchange_Value
          swap[1,j] = swap[1,j] + 1
        }
      }
      else if (runif(1)<=0.5){
        fraction_1 = log(g_tempered(Samples[[j+1]][i,1],1))-log((g_tempered(Samples[[j]][i,1],1)))
        alpha = min(1,exp(fraction_1*(-1/Temps[j+1]+1/Temps[j]))*
                      q_exchange(j+1,j,m)/q_exchange(j,j+1,m))
        if (runif(1)<alpha*c){
          Exchange_Value = Samples[[j]][i,1]
          Samples[[j]][i,1]=Samples[[j+1]][i,1]
          Samples[[j+1]][i,1]=Exchange_Value
          swap[1,j] = swap[1,j] + 1
        }
      }
      else {
        fraction_1 = log(g_tempered(Samples[[j-1]][i,1],1))-log((g_tempered(Samples[[j]][i,1],1)))
        alpha = min(1,exp(fraction_1*(-1/Temps[j-1]+1/Temps[j]))*
                      q_exchange(j-1,j,m)/q_exchange(j,j-1,m))
        if (runif(1)<alpha*c){
          Exchange_Value = Samples[[j]][i,1]
          Samples[[j]][i,1]=Samples[[j-1]][i,1]
          Samples[[j-1]][i,1]=Exchange_Value
          swap[1,j] = swap[1,j] + 1
        }
      }}}
  #gives percentage of swaps
  print (swap/n)
  return(Samples)}


set.seed(1234)
Samples = MH_Parallel_Tempering(sigma = c(0.001, 0.5, 0.6, 1, 5, 7, 10),x0 = c(30, 6, 27, 8, 15, 32, 28), n=1e4,
            temperatures = c(1, 2, 5, 10, 30, 50, 150))


#Visualisation
x_split = seq(0,45,by=0.01)
#Finding normalising constants
k = 1/(integrate(function(i)g(i),-Inf,3)$value + integrate(function(i)g(i),3,30)$value +
         integrate(function(i)g(i),30,Inf)$value )
k2 = 1/integrate(function(i)g_tempered(i,2),-Inf,Inf)$value
k3 = 1/integrate(function(i)g_tempered(i,5),-Inf,Inf)$value
k4 = 1/integrate(function(i)g_tempered(i,10),-Inf,Inf)$value
k5 = 1/integrate(function(i)g_tempered(i,30),-Inf,Inf)$value
k6 = 1/integrate(function(i)g_tempered(i,50),-Inf,Inf)$value
k7 = 1/integrate(function(i)g_tempered(i,150),-Inf,Inf)$value



plot(x_split,k*g_tempered(x_split,1),type="l",col="black",ylim=c(0,0.6),lty=1,ylab=
       "Density",xlab="x")
lines(x_split,k2*g_tempered(x_split,2),col="grey",lty=1,type="l")
lines(x_split,k3*g_tempered(x_split,5),col="yellow",lty=1,type="l")
lines(x_split,k4*g_tempered(x_split,10),col="blue",lty=1,type="l")
lines(x_split,k5*g_tempered(x_split,30),col="purple",lty=1,type="l")
lines(x_split,k6*g_tempered(x_split,50),col="orange",lty=1,type="l")
lines(x_split,k7*g_tempered(x_split,150),col="red",lty=1,type="l")
lines(density(Samples[[1]],1, adjust = 0.1), col = "green")
legend("topleft",legend=c("1","2","5","10","30","50","150"),col=c("black",
                                  "red","pink","orange","purple","blue","green"),lty=1,cex =0.75)



