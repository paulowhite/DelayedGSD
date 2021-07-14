#' @title Calculate boundaries for a group sequential design with delayed endpoints
#' @description Calculate boundaries (interim and decision) for a group sequential design with delayed endpoints based on planned and/or observed information using an error spending approach
#' 
#' @param kMax max number of analyses (including final)
#' @param sided one or two-sided
#' @param alpha type I error
#' @param beta type II error
#' @param InfoR.i planned or observed information rates at interim analysis, including the final analysis.
#' @param gammaA rho parameter for alpha error spending function
#' @param gammaB rho parameter for beta error spending function
#' @param method use method 1 or 2 from paper H&J
#' @param cNotBelowFixedc whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
#' @param delta effect that the study is powered for
#' @param InfoR.d (expected) information rate at each decision analysis (i.e. when stopping at an interim analysis). Should not include the final analysis.
#' @param bindingFutility whether the futility stopping rule is binding
#' @param trace whether to print some messages
#'
#' @examples
#' CalcBoundaries(kMax=2,
#'               sided=1,
#'               alpha=0.025,  
#'               beta=0.2,  
#'               InfoR.i=c(0.5,1),
#'               gammaA=2,
#'               gammaB=2,
#'               method=1,
#'               cNotBelowFixedc=TRUE,
#'               delta=1.5,
#'               InfoR.d=0.55)  
#' 
#' @export
CalcBoundaries <- function(kMax=2,  
                           sided=1,  
                           alpha=0.025, 
                           beta=0.2,  
                           InfoR.i=c(0.5,1),  
                           gammaA=2,  
                           gammaB=2,  
                           method=1, 
                           cNotBelowFixedc=FALSE, 
                           delta=1.5, 
                           InfoR.d=0.55,   
                           bindingFutility=TRUE,  #
                           trace=TRUE){  

    require(gsDesign)

    ## ** normalize user input
    call <- match.call() ## keep track of how the user run the function
    if(sided!=1){
        stop("Function cannot handle two-sided tests yet")
    }
  
    if(sum(InfoR.d < InfoR.i[-length(InfoR.i)])>0){
        stop("Information at decision analysis should not be smaller than information at interim analysis")
    }
  
    if(method %in% 1:2 == FALSE){
        stop("Please specify method=1 or method=2")
    }
    
    if(trace){message("In CalcBoundaries, the method assumes that positive effects are good")}
 
    ## ** compute boundaries at interim
    if(bindingFutility){
        StandardDesign <- gsDesign::gsDesign(k=kMax, test.type=3,alpha=alpha,beta=beta,
                                             timing=InfoR.i,
                                             sfu=sfPower,sfupar=gammaA,
                                             sfl=sfPower,sflpar=gammaB)
    } else {
        StandardDesign <- gsDesign::gsDesign(k=kMax, test.type=4,alpha=alpha,beta=beta,
                                             timing=InfoR.i,
                                             sfu=sfPower,sfupar=gammaA,
                                             sfl=sfPower,sflpar=gammaB)
    }
  
  
                                        #R <- getDesignCharacteristics(StandardDesign)$inflationFactor
    R <- StandardDesign$n.I[kMax]
    Info.max <- ((qnorm(1-alpha)+qnorm(1-beta))/delta)^2*R
    Info.d <- InfoR.d*Info.max
  
    if(max(Info.d) >= Info.max){
        stop("Function cannot handle Id >= Imax yet")
    }
    
                                        #uk <- StandardDesign$criticalValues
                                        #lk <- c(StandardDesign$futilityBounds,uk[kMax])
    uk <- StandardDesign$upper$bound
  
  #browser()

    ## ** compute boundaries at decision and possibly update futility boundary at interim
    lk <- StandardDesign$lower$bound
    ck <- rep(0,kMax-1)
    cMin <- ifelse(cNotBelowFixedc,qnorm(1-alpha),-Inf)

    if(method==1){
    
        for(k in 1:(kMax-1)){
            ck[k] <- Method1(uk = uk[1:k],
                             lk = lk[1:k],
                             Ik = (InfoR.i*Info.max)[1:k],
                             Id = Info.d[k],
                             Imax = Info.max,
                             cMin = cMin,
                             bindingFutility = bindingFutility)
        }
    
    } else if(method==2){
    
        for(k in 1:(kMax-1)){
            delayedBnds <- Method2(uk = uk[1:k],
                                   lk = lk[1:k],
                                   Ik = (InfoR.i*Info.max)[1:k],
                                   Id = Info.d[k],
                                   Imax = Info.max,
                                   delta = delta,
                                   cMin = cMin,
                                   bindingFutility = bindingFutility)
            ck[k] <- delayedBnds$critval
            lk[k] <- delayedBnds$lkd
        }
    
    }


    ## ** output
    out <- list(call = call,
                uk = uk, 
                lk = lk,
                ck = ck,
                Info.i = InfoR.i*Info.max,
                Info.d = Info.d,
                Info.max = Info.max,
                InflationFactor = R,
                alpha = alpha,
                kMax = kMax,
                sided = sided,
                beta = beta,
                gammaA = gammaA,
                gammaB = gammaB,
                method = method,
                delta = delta,
                cNotBelowFixedc = cNotBelowFixedc,
                cMin = cMin,
                bindingFutility = bindingFutility)
    class(out) <- append("delayedGSD",class(out))
    return(out)
}

print.delayedGSD <- function(x, ...){
  
  print(x$call)
  
  cat("\n")
  
  bnds <- cbind(c(1:length(x$lk)),x$lk,x$uk,c(x$ck,x$uk[length(x$uk)]))
  colnames(bnds) <- c("Stage","Lower","Upper","Decision")
  
  cat("Boundaries \n")
  print(bnds)
  
  cat("\n")
  
  cat("Planned maximum information \n")
  print(x$Info.max)
  
  cat("\n")
  
  cat("Inflation factor \n")
  print(x$InflationFactor)
  
  
}
