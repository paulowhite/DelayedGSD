#calculate critval c using method 1:
require(BB)

Method1 <- function(uk,  #upper bounds for all analyses up to and including current stage k
                    lk,  #lower bounds for all analyses up to and including current stage k
                    Ik,  #Information for all analyses up to and including current stage k
                    Id,  #Observed information at decision analysis k
                    Imax,  #Maximum information
                    sided=1,  #one or two sided
                    cMin=-Inf, # minimun possible value c for the decision analysis, typically that for a fixed sample test (H & J page 10)
                    bindingFutility=TRUE){  #whether the futility stopping rule is binding
  
    require(mvtnorm)
    
    if(sided!=1){
        stop("Function cannot handle two-sided tests yet")
    }
  
    if(Id >= Imax){
        stop("Function cannot handle Id >= Imax yet")
    }
  
    ## message("the method assumes that positive effects are good")
  
    k <- length(uk)
    Ik <- c(Ik,Id)
    sigmaZk <- diag(1,k+1)
    for(i in 1:(k+1)){
        for(j in i:(k+1)){
            sigmaZk[i,j] <- sqrt(Ik[i]/Ik[j])
            sigmaZk[j,i] <- sqrt(Ik[i]/Ik[j])
        }
    }
  
    if(bindingFutility){
      f <- function(x){
          y <- pmvnorm(lower = c(lk[0:(k-1)],uk[k],-Inf),
                       upper = c(uk[0:(k-1)],Inf,x),
                       mean=rep(0,k+1),
                       sigma= sigmaZk) - 
              pmvnorm(lower = c(lk[0:(k-1)],-Inf,x),
                      upper = c(uk[0:(k-1)],lk[k],Inf),
                      mean=rep(0,k+1),
                      sigma= sigmaZk)
          as.numeric(y)
      }
    } else {
      lk2 <- rep(-Inf,length(lk))
      
      f <- function(x){
        y <- pmvnorm(lower = c(lk2[0:(k-1)],uk[k],-Inf),
                     upper = c(uk[0:(k-1)],Inf,x),
                     mean=rep(0,k+1),
                     sigma= sigmaZk) - 
          pmvnorm(lower = c(lk2[0:(k-1)],-Inf,x),
                  upper = c(uk[0:(k-1)],lk[k],Inf),
                  mean=rep(0,k+1),
                  sigma= sigmaZk)
        as.numeric(y)
      }
    }
  
    c <- try(uniroot(f, lower = lk[k], upper = uk[k]*1.1)$root,silent=T)
    if(inherits(c,"try-error")){
        warning("uniroot produced an error, switching to dfsane")
        c <- BB::dfsane(f, par = (uk[k]*1.1+lk[k])/2, control = list(maxit = 1000, tol = 1e-7), quiet = TRUE)$par
    }
    max(c,cMin)
}
