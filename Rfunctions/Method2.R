#calculate critval c using method 2:

Method2 <- function(uk,  #upper bounds for all analyses up to and including current stage k
                    lk,  #lower bounds for all analyses up to and including current stage k
                    Ik,  #Information for all analyses up to and including current stage k
                    Id,  #Expected information at decision analysis
                    Imax, #planned maximum information
                    delta, #planned effect
                    sided=1,  #whether the test is 1 or 2-sided
                    cMin=-Inf){  # minimun possible value c for the decision analysis, typically that for a fixed sample test (H & J page 10)
  
  require(mvtnorm)
  
  if(sided!=1){
    stop("Function cannot handle two-sided tests yet")
  }
  
  ## message("the method assumes that positive effects are good")
  
  if(Id >= Imax){
    stop("Function cannot handle Id >= Imax yet")
  }
  
  k <- length(uk)
  Ik <- c(Ik,Id)
  sigmaZk <- diag(1,k+1)
  for(i in 1:(k+1)){
    for(j in i:(k+1)){
      sigmaZk[i,j] <- sqrt(Ik[i]/Ik[j])
      sigmaZk[j,i] <- sqrt(Ik[i]/Ik[j])
    }
  }
  
  theta <- delta*sqrt(Ik)
  
  #browser()
  
  sol <- uniroot(function(x){
    
    c <- Method1(lk=c(lk[0:(k-1)],x),uk=uk,Ik=Ik[1:k],Id=Id,Imax=Imax,cMin=cMin)
    
    pmvnorm(lower = c(lk[0:(k-1)],uk[k],-Inf),
            upper = c(uk[0:(k-1)],Inf,c),
            mean=theta,
            sigma= sigmaZk) + 
      pmvnorm(lower = c(lk[0:(k-1)],-Inf,-Inf),
              upper = c(uk[0:(k-1)],x,c),
              mean=theta,
              sigma= sigmaZk) - (ErrorSpend(I=Id,Imax=Imax,beta_or_alpha=0.2,rho=2)-ErrorSpend(I=Ik[k],Imax=Imax,beta_or_alpha=0.2,rho=2))},
    lower = -10,
    upper = uk[k]*0.999)$root #can't be exactly uk[k] since upper bound for method1 function and then the lower and upper bound would be equal
  
  c <- Method1(lk=c(lk[0:(k-1)],sol),uk=uk,Ik=Ik[1:k],Id=Id,Imax=Imax,cMin=cMin)
  
  list(lkd=sol,critval=c)
}
