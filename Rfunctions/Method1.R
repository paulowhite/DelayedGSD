#calculate critval c using method 1:

Method1 <- function(uk,  #upper bounds for all analyses up to and including current stage k
                    lk,  #lower bounds for all analyses up to and including current stage k
                    Ik,  #Information for all analyses up to and including current stage k
                    Id,  #Observed information at decision analysis k
                    Imax,  #Maximum information
                    sided=1){ #one or two sided
  
  require(mvtnorm)
  
  if(sided!=1){
    stop("Function cannot handle two-sided tests yet")
  }
  
  if(Id >= Imax){
    stop("Function cannot handle Id >= Imax yet")
  }
  
  message("the method assumes that positive effects are good")
  
  k <- length(uk)
  Ik <- c(Ik,Id)
  sigmaZk <- diag(1,k+1)
  for(i in 1:(k+1)){
    for(j in i:(k+1)){
      sigmaZk[i,j] <- sqrt(Ik[i]/Ik[j])
      sigmaZk[j,i] <- sqrt(Ik[i]/Ik[j])
    }
  }
  
  c <- uniroot(function(x){pmvnorm(lower = c(lk[0:(k-1)],uk[k],-Inf),
                                   upper = c(uk[0:(k-1)],Inf,x),
                                   mean=rep(0,k+1),
                                   sigma= sigmaZk) - 
      pmvnorm(lower = c(lk[0:(k-1)],-Inf,x),
              upper = c(uk[0:(k-1)],lk[k],Inf),
              mean=rep(0,k+1),
              sigma= sigmaZk)},
      lower = lk[k],
      upper = uk[k]*1.1)$root
  c
}
