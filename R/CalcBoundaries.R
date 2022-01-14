## * CalcBoundaries (documentation)
#' @title Calculate boundaries for a group sequential design with delayed endpoints
#' @description Calculate boundaries (interim and decision) for a group sequential design with delayed endpoints based on planned and/or observed information using an error spending approach
#' 
#' @param kMax max number of analyses (including final)
#' @param sided one or two-sided
#' @param alpha type I error
#' @param beta type II error
#' @param InfoR.i planned or observed information rates at interim analysis, including the final analysis.
#' @param rho_alpha rho parameter for alpha error spending function
#' @param rho_beta rho parameter for beta error spending function
#' @param method use method 1 or 2 from paper H&J
#' @param cNotBelowFixedc whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
#' @param delta effect that the study is powered for
#' @param InfoR.d (expected) information rate at each decision analysis (i.e. when stopping at an interim analysis). Should not include the final analysis.
#' @param bindingFutility whether the futility stopping rule is binding
#' @param alternative a character string specifying the alternative hypothesis, "greater" or "less".
#' greater is for Ho= theta > 0, "less" is for Ho= theta < 0 (note that in Jennison and Turnbull's book chapter (2013) they consider less)
#' @param n planned sample size in each group. Optional argument.
#' @param trace whether to print some messages
#'
#' @examples
#' myBound <- CalcBoundaries(kMax=2,
#'               sided=1,
#'               alpha=0.025,  
#'               beta=0.2,  
#'               InfoR.i=c(0.5,1),
#'               rho_alpha=2,
#'               rho_beta=2,
#'               method=2,
#'               cNotBelowFixedc=TRUE,
#'               delta=1.5,
#'               InfoR.d=0.55)
#' myBound
#' plot(myBound)
#' 
#' ## to reproduce bounds from CJ DSBS course slide 106
#' myBound <- CalcBoundaries(kMax=3,
#'               sided=1,
#'               alpha=0.025,  
#'               beta=0.1,  
#'               InfoR.i=c(3.5,6.75,12)/12,
#'               rho_alpha=1.345,
#'               rho_beta=1.345,
#'               method=1, ## has been changed from 2 to 1
#'               cNotBelowFixedc=TRUE,
#'               bindingFutility=FALSE,
#'               delta=1,
#'               InfoR.d=c(5.5,8.75)/12)
#' myBound
#' plot(myBound)

## * CalcBoundaries (code)
#' @export
CalcBoundaries <- function(kMax=2,  
                           sided=1,  
                           alpha=0.025, 
                           beta=0.2,  
                           InfoR.i=c(0.5,1),  
                           rho_alpha=2,  
                           rho_beta=2,  
                           method=1, 
                           cNotBelowFixedc=FALSE, 
                           delta=1.5, 
                           InfoR.d=0.55,   
                           bindingFutility=TRUE,
                           alternative = "less",
                           n=NULL,
                           trace=TRUE){  

    requireNamespace("gsDesign")

    ## ** normalize user input
    call <- match.call() ## keep track of how the user run the function
    alternative <- match.arg(alternative, c("less","greater"))
    if(sided!=1){
        stop("Function cannot handle two-sided tests yet")
    }
  
    if(sum(InfoR.d < InfoR.i[-length(InfoR.i)])>0){
        stop("Information at decision analysis should not be smaller than information at interim analysis")
    }
  
    if(method %in% 1:3 == FALSE){
        stop("Please specify method=1, method=2, or method=3.")
    }
    
    if( (alternative=="less" & delta<0) | (alternative=="greater" & delta>0)){
        stop("The values given for arguments alternative and delta are inconsistent.\n When alternative=less, delta should be positive.\n When alternative=greater, delta should be negative.")
    }

    ## ** compute boundaries at decision and possibly update futility boundary at interim
    cMin <- ifelse(cNotBelowFixedc,stats::qnorm(1-alpha),-Inf)
    
    ## ## ** remove boundaries corresponding to stage that will not be reached
    ## ## e.g. we stop early (stage 1) and need to re-compute the boundary at decision
    ## ##      futility and efficacy boundaries at stage 2 or later will never be used
    ## ##      same for the decision boundaries at stage 2 or later

    ## indexNNA.i <- which(!is.na(InfoR.i))
    ## indexNNA.d <- which(!is.na(InfoR.d))
    ## uk[-indexNNA.i] <- NA
    ## lk[-indexNNA.i] <- NA
    ## ck[-indexNNA.d] <- NA

    ## if(any(max(Info.d[indexNNA.d]) >= Info.max)){
    ##     warning("Information at decision exceed maximum planned information. \n",
    ##             "Use the maximum planned information to compute the boundary at decision. \n")
    ##     browser()
    ##     ## max(Info.d[indexNNA.d]) >= Info.max)
    ## }


    if(method==1){
    
    if(bindingFutility){test.type <- 3}else{test.type <- 4}

        if(alternative=="greater"){
            delta <- -delta
        }
        ## ** compute boundaries at interim
        StandardDesign <- gsDesign::gsDesign(k=kMax, test.type=test.type, alpha=alpha, beta=beta,
                                             timing=InfoR.i,
                                             sfu=gsDesign::sfPower, sfupar=rho_alpha,
                                             sfl=gsDesign::sfPower, sflpar=rho_beta)
        ## R <- getDesignCharacteristics(StandardDesign)$inflationFactor
        R <- StandardDesign$n.I[kMax]
        Info.max <- ((stats::qnorm(1-alpha)+stats::qnorm(1-beta))/delta)^2*R
        Info.d <- InfoR.d*Info.max

        ##uk <- StandardDesign$criticalValues
        ##lk <- c(StandardDesign$futilityBounds,uk[kMax])
        uk <- StandardDesign$upper$bound
        lk <- StandardDesign$lower$bound
        ck <- rep(0,kMax-1)
        alphaSpent <- StandardDesign$upper$spend
        betaSpent <- StandardDesign$lower$spend
        for(k in 1:(kMax-1)){
            ck[k] <- Method1(uk = uk[1:k],
                             lk = lk[1:k],
                             Info.i = (InfoR.i*Info.max)[1:k],
                             Info.d = Info.d[k],
                             Info.max = Info.max,
                             cMin = cMin,
                             bindingFutility = bindingFutility)
        }

        if(alternative=="greater"){
            uk <- -uk
            lk <- -lk
            ck <- -ck
            delta <- -delta
        }
        
    } else if(method==2){
    
        delayedBnds <- Method2(rho_alpha = rho_alpha,
                               rho_beta = rho_beta,
                               alpha = alpha,
                               beta = beta, 
                               Kmax = kMax,
                               Info.max = NULL,
                               InfoR.i = InfoR.i,
                               InfoR.d = InfoR.d,
                               delta = delta, 
                               alternative = alternative,
                               Trace = trace,
                               cMin = cMin)
        uk <- delayedBnds$boundaries$u.k
        lk <- delayedBnds$boundaries$l.k
        ck <- delayedBnds$boundaries$c.k[1:(kMax-1)]
        R <- delayedBnds$coef
        Info.d <- delayedBnds$Info.d
        Info.max <- delayedBnds$Info.max
        alphaSpent <- delayedBnds$boundaries$Inc.Type.I
        betaSpent <- delayedBnds$boundaries$Inc.Type.II

    }else if(method==3){

        delayedBnds <- Method3(rho_alpha = rho_alpha,
                               rho_beta = rho_beta,
                               alpha = alpha,
                               beta = beta, 
                               Kmax = kMax,
                               Info.max = NULL,
                               InfoR.i = InfoR.i,
                               InfoR.d = InfoR.d,
                               delta = delta, 
                               alternative = alternative,
                               Trace = trace) 

        uk <- delayedBnds$boundaries$u.k
        lk <- delayedBnds$boundaries$l.k
        ck <- delayedBnds$boundaries$c.k[1:(kMax-1)]
        R <- delayedBnds$coef
        Info.d <- delayedBnds$Info.d
        Info.max <- delayedBnds$Info.max
        alphaSpent <- delayedBnds$boundaries$Inc.Type.I
        betaSpent <- delayedBnds$boundaries$Inc.Type.II

    }

    ## ** output
    out <- list(call = call,
                n.obs = n,
                stage = data.frame(k = 0, type = "planning"),
                conclusion = matrix(as.character(NA), nrow = 4, ncol = kMax, dimnames = list(c("interim","reason.interim","decision","comment.decision"),NULL)),
                uk = uk, 
                lk = lk,
                ck = ck,
                Info.i = InfoR.i*Info.max,
                Info.d = Info.d,
                Info.max = Info.max,
                lmm = vector(mode = "list", length = kMax),
                InflationFactor = R,
                alpha = alpha,
                alphaSpent = cumsum(alphaSpent),
                kMax = kMax,
                sided = sided,
                beta = beta,
                betaSpent = cumsum(betaSpent),
                rho_alpha = rho_alpha,
                rho_beta = rho_beta,
                method = method,
                delta = delta,
                cNotBelowFixedc = cNotBelowFixedc,
                cMin = cMin,
                bindingFutility = bindingFutility,
                alternative = alternative,
                planned = list(uk = uk,
                               lk = lk,
                               ck = ck,
                               alphaSpent = cumsum(alphaSpent),
                               betaSpent = cumsum(betaSpent),
                               Info.i = InfoR.i*Info.max,
                               Info.d = Info.d,
                               delta = delta))
    class(out) <- append("delayedGSD",class(out))
    return(out)
}


