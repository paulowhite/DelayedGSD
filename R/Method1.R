#' @title calculate critval c using method 1
#'
#' @param uk upper bounds for all analyses up to and including current stage k
#' @param lk lower bounds for all analyses up to and including current stage k
#' @param Info.i Information for all analyses up to and including current stage k
#' @param Info.d Observed information at decision analysis k
#' @param Info.max Maximum information
#' @param cMin minimun possible value c for the decision analysis, typically that for a fixed sample test (H & J page 10)
#' @param ImaxAnticipated set to TRUE if c should be calculated according to Eq 15 in HJ because study was stopped early due to Imax reached
#' @param rho_alpha value of rho used for type I error spending function
#' @param alpha one-sided alpha level to be used for the study
#' @param bindingFutility whether the futility stopping rule is binding
#'
#' @export
Method1 <- function(uk,  
                    lk,  
                    Info.i,
                    Info.d,
                    Info.max=NULL,
                    cMin=-Inf,
                    ImaxAnticipated=FALSE,
                    rho_alpha=2, 
                    alpha=0.025, 
                    bindingFutility=TRUE){  
  
    requireNamespace("mvtnorm")
    requireNamespace("BB")
  
    ## message("the method assumes that positive effects are good")
  
    k <- length(uk)
    
    ## If interim analysis k was skipped due to max info reached at interim or expected to be reached at decision, then calculate c according to Eq 15 from HJ
    if(ImaxAnticipated==TRUE){
                                        #stop("Function cannot handle Info.d >= Imax yet")
      
      Info.all <- c(Info.i[1:(k-1)],Info.d)
      sigmaZk <- diag(1,k)
      for(i in 1:(k)){
        for(j in i:(k)){
          sigmaZk[i,j] <- sqrt(Info.all[i]/Info.all[j])
          sigmaZk[j,i] <- sqrt(Info.all[i]/Info.all[j])
        }
      }
      
        if(bindingFutility){
            f <- function(x){
                y <- mvtnorm::pmvnorm(lower = c(lk[0:(k-1)],x),
                                      upper = c(uk[0:(k-1)],Inf),
                                      mean=rep(0,k),
                                      sigma= sigmaZk) - (alpha-ErrorSpend(I = Info.all[k-1], rho = rho_alpha,beta_or_alpha = alpha, Info.max = Info.max))
                as.numeric(y)
            }
        } else {
            lk2 <- rep(-Inf,length(lk))
        
            f <- function(x){
                y <- mvtnorm::pmvnorm(lower = c(lk2[0:(k-1)],x),
                                      upper = c(uk[0:(k-1)],Inf),
                                      mean=rep(0,k),
                                      sigma= sigmaZk) - (alpha-ErrorSpend(I = Info.all[k-1], rho = rho_alpha,beta_or_alpha = alpha, Info.max = Info.max))
                as.numeric(y)
            }
        }
        ## usual way of computing c based on Eq 14 HJ  
    } else {
      
      Info.all <- c(Info.i,Info.d)
      sigmaZk <- diag(1,k+1)
      for(i in 1:(k+1)){
        for(j in i:(k+1)){
          sigmaZk[i,j] <- sqrt(Info.all[i]/Info.all[j])
          sigmaZk[j,i] <- sqrt(Info.all[i]/Info.all[j])
        }
      }
      
        if(bindingFutility){
            f <- function(x){
                y <- mvtnorm::pmvnorm(lower = c(lk[0:(k-1)],uk[k],-Inf),
                                      upper = c(uk[0:(k-1)],Inf,x),
                                      mean=rep(0,k+1),
                                      sigma= sigmaZk) - 
                    mvtnorm::pmvnorm(lower = c(lk[0:(k-1)],-Inf,x),
                                     upper = c(uk[0:(k-1)],lk[k],Inf),
                                     mean=rep(0,k+1),
                                     sigma= sigmaZk)
                as.numeric(y)
            }
        } else {
            lk2 <- rep(-Inf,length(lk))
        
            f <- function(x){
                y <- mvtnorm::pmvnorm(lower = c(lk2[0:(k-1)],uk[k],-Inf),
                                      upper = c(uk[0:(k-1)],Inf,x),
                                      mean=rep(0,k+1),
                                      sigma= sigmaZk) - 
                    mvtnorm::pmvnorm(lower = c(lk2[0:(k-1)],-Inf,x),
                                     upper = c(uk[0:(k-1)],lk[k],Inf),
                                     mean=rep(0,k+1),
                                     sigma= sigmaZk)
                as.numeric(y)
            }
        }
    }

    c <- try(stats::uniroot(f,
                            lower = lk[utils::tail(intersect(which(!is.infinite(lk)),1:k),1)], ## last boundary among the k already computed that is not infinite  
                            upper = uk[utils::tail(intersect(which(!is.infinite(uk)),1:k),1)]*1.1)$root,
             silent=T)
    if(inherits(c,"try-error")){
        warning("uniroot produced an error, switching to dfsane")
        c <- BB::dfsane(f,
                        par = (uk[utils::tail(intersect(which(!is.infinite(uk)),1:k),1)]*1.1+lk[utils::tail(intersect(which(!is.infinite(lk)),1:k),1)])/2,
                        control = list(maxit = 1000, tol = 1e-7), quiet = TRUE)$par
    }
    max(c,cMin)
}
