#' @title Update c_k based on observed information at the decision analysis for Method 2
#' @description Update c_k based on observed information at the decision analysis. The value will be the minimum of the re-calculated value for c_k based on the observed information and cMin as from the boundary object that is provided as input.
#' 
#' 
#' 
#' @param obj       output from either NonBindingHJ.R or Method2_PC.R
#' @param k         the analysis k at which to recalculate c_k  
#' @param Info.new  the observed information at decision analysis k
#'
#' @examples
#' bCJ <- NonBindingHJ(rho_alpha=1.345,
#'                     rho_beta=1.345,
#'                     alpha=0.025,
#'                     beta=0.1,
#'                     Kmax=3,
#'                     Info.max=12,
#'                     InfoR.i=c(3.5,6.75,12)/12,
#'                     InfoR.d=c(5.5,8.75)/12,
#'                     delta=1,  
#'                     abseps = 1e-06, 
#'                     direction="smaller",
#'                     sided=1
#'                    )
#'
#' Update_ck(bCJ,k=2,Info.new=8)

Update_ck <- function(obj,       #Output from NonBindingHJ.R
                      k,         #the analysis at which to update c
                      Info.new){ #the new observed information at decision analysis k
  
  Info.i <- obj$Info.i
  Kmax <- obj$Kmax
  lk <- obj$boundaries[,"l.k"]
  uk <- obj$boundaries[,"u.k"]
  IncAlpha <- obj$boundaries[,"Inc.Type.I"]
  binding <- obj$binding
  
  if(obj$direction=="greater"){
    lk <- -lk
    uk <- -uk
  }
  
  sigmaZk <- diag(1,Kmax)
  for(i in 1:Kmax){
    for(j in i:Kmax){
      sigmaZk[i,j] <- sqrt(Info.i[i]/Info.i[j])
      sigmaZk[j,i] <- sqrt(Info.i[i]/Info.i[j])
    }
  }
  
  sigmaZk2 <- matrix(NA,ncol=k+1,nrow=k+1)
  sigmaZk2[1:k,1:k] <- sigmaZk[1:k,1:k]
  sigmaZk2[k+1,k+1] <- 1
  sigmaZk2[1:k,k+1] <- sigmaZk2[k+1,1:k] <- sqrt(Info.i[1:k]/Info.new)
  
  if(binding){
    stop("this has not been implemented yet")
  } else {
    c_new <- uniroot(function(x){pmvnorm(lower = c(rep(-Inf,k-1),uk[k],x),
                                         upper = c(uk[1:(k-1)],Inf,Inf),
                                         mean=rep(0,k+1),
                                         sigma= sigmaZk2,
                                         abseps = obj$abseps) - IncAlpha[k]},
                     lower = lk[k-1],
                     upper = uk[k-1],
                     tol = obj$abseps)$root
  }
  
  max(obj$cMin,c_new)
}

