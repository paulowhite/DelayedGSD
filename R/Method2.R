#' @title calculate critval c using method 2
#'
#' @param uk upper bounds for all analyses up to and including current stage k
#' @param lk lower bounds for all analyses up to and including current stage k
#' @param Info.i Information for all analyses up to and including current stage k
#' @param Info.d Observed information at decision analysis k
#' @param Info.max Maximum information
#' @param delta planned effect
#' @param sided one or two sided
#' @param cMin minimun possible value c for the decision analysis, typically that for a fixed sample test (H & J page 10)
#' @param bindingFutility whether the futility stopping rule is binding
#'
#' @export
Method2 <- function(uk, 
                    lk, 
                    Info.i,
                    Info.d, 
                    Info.max,
                    delta,
                    sided=1,
                    cMin=-Inf, 
                    bindingFutility=TRUE){ 
  
  require(mvtnorm)
  
  if(sided!=1){
    stop("Function cannot handle two-sided tests yet")
  }
  
  if(!bindingFutility){
    stop("Function cannot handle non-binding futility for Method 2")
  }
  
  ## message("the method assumes that positive effects are good")
  
  if(Info.d >= Info.max){
    stop("Function cannot handle Info.d >= Info.max yet")  ##UPDATE, SEE METHOD 1 WHERE THIS IS HANDLED
  }
  
  k <- length(uk)
  Info.all <- c(Ik,Id)
  sigmaZk <- diag(1,k+1)
  for(i in 1:(k+1)){
    for(j in i:(k+1)){
      sigmaZk[i,j] <- sqrt(Info.all[i]/Info.all[j])
      sigmaZk[j,i] <- sqrt(Info.all[i]/Info.all[j])
    }
  }
  
  theta <- delta*sqrt(Info.all)
  
  #browser()
  
  sol <- uniroot(function(x){
    
    c <- Method1(lk=c(lk[0:(k-1)],x),uk=uk,Info.i=Info.all[1:k],Info.d=Info.d,Info.max=Info.max,cMin=cMin)
    
    pmvnorm(lower = c(lk[0:(k-1)],uk[k],-Inf),
            upper = c(uk[0:(k-1)],Inf,c),
            mean=theta,
            sigma= sigmaZk) + 
      pmvnorm(lower = c(lk[0:(k-1)],-Inf,-Inf),
              upper = c(uk[0:(k-1)],x,c),
              mean=theta,
              sigma= sigmaZk) - (ErrorSpend(I=Info.d,Imax=Info.max,beta_or_alpha=0.2,rho=2)-ErrorSpend(I=Info.all[k],Imax=Info.max,beta_or_alpha=0.2,rho=2))},
    lower = -10,
    upper = uk[k]*0.999)$root #can't be exactly uk[k] since upper bound for method1 function and then the lower and upper bound would be equal
  
  c <- Method1(lk=c(lk[0:(k-1)],sol),uk=uk,Info.i=Info.i[1:k],Info.d=Info.d,Info.max=Info.max,cMin=cMin)
  
  list(lkd=sol,critval=c)
}
