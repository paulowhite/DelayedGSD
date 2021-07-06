##' @title add title
##' @description add description
##' 
##' @param kMax max number of analyses (including final)
##' @param sided one or two-sided
##' @param alpha type I error
##' @param beta type II error
##' @param InfoR.i planned or observed information rates
##' @param gammaA rho parameter for alpha error spending function
##' @param gammaB rho parameter for beta error spending function
##' @param method use method 1 or 2 from paper H&J
##' @param cNotBelowFixedc whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
##' @param delta effect that the study is powered for
##' @param InfoR.d (expected) information rate at each decision analysis
##' @param bindingFutility whether the futility stopping rule is binding
##' @param Trace whether to print some messages
##'
##' @examples
##' ## add example
##' 
##' @export
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
                           Trace=TRUE){  
                            
  require(rpact)

    ## ** normalize user input
    if(sided!=1){
        stop("Function cannot handle two-sided tests yet")
    }
  
    if(sum(InfoR.d < InfoR.i[-length(InfoR.i)])>0){
        stop("Information at decision analysis should not be smaller than information at interim analysis")
    }
  
    if(method %in% 1:2 == FALSE){
        stop("Please specify method=1 or method=2")
    }
    
    if(Trace){message("In CalcBoundaries, the method assumes that positive effects are good")}
 
    ## ** compute boundaries at interim
    if(bindingFutility){
        StandardDesign <- gsDesign(k=kMax, test.type=3,alpha=alpha,beta=beta,
                                   timing=InfoR.i,
                                   sfu=sfPower,sfupar=gammaA,
                                   sfl=sfPower,sflpar=gammaB)
    } else {
        StandardDesign <- gsDesign(k=kMax, test.type=4,alpha=alpha,beta=beta,
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

    if(method==1){
    
        for(k in 1:(kMax-1)){
            ck[k] <- Method1(uk = uk[1:k],
                             lk = lk[1:k],
                             Ik = (InfoR.i*Info.max)[1:k],
                             Id = Info.d[k],
                             Imax = Info.max,
                             cMin = ifelse(cNotBelowFixedc,qnorm(1-alpha),-Inf),
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
                                   cMin = ifelse(cNotBelowFixedc,qnorm(1-alpha),-Inf),
                                   bindingFutility = bindingFutility)
            ck[k] <- delayedBnds$critval
            lk[k] <- delayedBnds$lkd
        }
    
    }


    ## ** output
    out <- list(uk = uk, 
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
                bindingFutility = bindingFutility)
    class(out) <- append("delayedGSD",class(out))
    return(out)
}

print.delayedGSD <- function(x, ...){
}
