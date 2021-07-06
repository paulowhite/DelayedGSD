#calculate critval c using method 1:
require(BB)

Method1 <- function(uk,  #upper bounds for all analyses up to and including current stage k
                    lk,  #lower bounds for all analyses up to and including current stage k
                    Info.i,  #Information for all analyses up to and including current stage k
                    Info.d,  #Observed information at decision analysis k
                    Info.max=NULL,  #Maximum information
                    sided=1,  #one or two sided
                    cMin=-Inf, # minimun possible value c for the decision analysis, typically that for a fixed sample test (H & J page 10)
                    ImaxAnticipated=FALSE, #set to TRUE if c should be calculated according to Eq 15 in HJ because study was stopped early due to Imax reached
                    rho=2, #value of rho used for type I error spending function
                    alpha=0.025, #one-sided alpha level to be used for the study
                    bindingFutility=TRUE){  #whether the futility stopping rule is binding
  
    require(mvtnorm)
    
    if(sided!=1){
        stop("Function cannot handle two-sided tests yet")
    }
  
    ## message("the method assumes that positive effects are good")
  
    k <- length(uk)
    
    #If interim analysis k was skipped due to max info reached at interim or expected to be reached at decision, then calculate c according to Eq 15 from HJ
    if(ImaxAnticipated=TRUE){
      #stop("Function cannot handle Info.d >= Imax yet")
      
      Info.i <- c(Info.i[1:(k-1)],Info.d)
      sigmaZk <- diag(1,k)
      for(i in 1:(k)){
        for(j in i:(k)){
          sigmaZk[i,j] <- sqrt(Info.i[i]/Info.i[j])
          sigmaZk[j,i] <- sqrt(Info.i[i]/Info.i[j])
        }
      }
      
      if(bindingFutility){
        f <- function(x){
          y <- pmvnorm(lower = c(lk[0:(k-1)],x),
                       upper = c(uk[0:(k-1)],Inf),
                       mean=rep(0,k),
                       sigma= sigmaZk) - (alpha-ErrorSpend(I = Info.i[k-1], rho = rho,beta_or_alpha = alpha, Info.max = Info.max))
          as.numeric(y)
        }
      } else {
        lk2 <- rep(-Inf,length(lk))
        
        f <- function(x){
          y <- pmvnorm(lower = c(lk2[0:(k-1)],x),
                       upper = c(uk[0:(k-1)],Inf),
                       mean=rep(0,k),
                       sigma= sigmaZk) - (alpha-ErrorSpend(I = Info.i[k-1], rho = rho,beta_or_alpha = alpha, Info.max = Info.max))
          as.numeric(y)
        }
      }
    #usual way of computing c based on Eq 14 HJ  
    } else {
      
      Info.i <- c(Info.i,Info.d)
      sigmaZk <- diag(1,k+1)
      for(i in 1:(k+1)){
        for(j in i:(k+1)){
          sigmaZk[i,j] <- sqrt(Info.i[i]/Info.i[j])
          sigmaZk[j,i] <- sqrt(Info.i[i]/Info.i[j])
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
    }
  
    c <- try(uniroot(f, lower = lk[k], upper = uk[k]*1.1)$root,silent=T)
    if(inherits(c,"try-error")){
        warning("uniroot produced an error, switching to dfsane")
        c <- BB::dfsane(f, par = (uk[k]*1.1+lk[k])/2, control = list(maxit = 1000, tol = 1e-7), quiet = TRUE)$par
    }
    max(c,cMin)
}
