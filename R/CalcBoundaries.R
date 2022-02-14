## * CalcBoundaries (documentation)
#' @title Calculate boundaries for a group sequential design with delayed endpoints
#' @description Calculate boundaries (interim and decision) for a group sequential design with delayed endpoints based on planned and/or observed information using an error spending approach
#' 
#' @param kMax max number of analyses (including final)
#' @param alpha type I error
#' @param beta type II error
#' @param InfoR.i planned or observed information rates at interim analysis, including the final analysis.
#' @param rho_alpha rho parameter for alpha error spending function
#' @param rho_beta rho parameter for beta error spending function
#' @param method use method 1 or 2 from paper H&J
#' @param cNotBelowFixedc whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
#' @param delta effect that the study is powered for
#' @param InfoR.d (expected) information rate at each decision analysis (i.e. when stopping at an interim analysis). Should not include the final analysis.
#' @param bindingFutility [logical]  whether the futility stopping rule is binding.
#' @param alternative a character string specifying the alternative hypothesis, \code{"greater"} or \code{"less"}.
#' H0 \eqn{\theta=0} vs H1 \eqn{theta<0} (\code{"less"}) or theta > 0 (\code{"greater"}).
#' Note that in Jennison and Turnbull's book chapter (2013) they consider greater.
#' @param n planned sample size in each group. Optional argument.
#' @param trace whether to print some messages
#'
#' @examples
#' myBound <- CalcBoundaries(kMax=2,
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
                           alternative = "greater",
                           n=NULL,
                           trace=FALSE){  

    requireNamespace("gsDesign")

    ## ** normalize user input
    call <- match.call() ## keep track of how the user run the function
    alternative <- match.arg(alternative, c("less","greater"))
  
    if(sum(InfoR.d < InfoR.i[-length(InfoR.i)])>0){
        stop("Information at decision analysis should not be smaller than information at interim analysis")
    }
  
    if(method %in% 1:3 == FALSE){
        stop("Please specify method=1, method=2, or method=3.")
    }

    if(alternative=="greater" & delta<0){
        stop("The values given for arguments \'alternative\' and \'delta\' are inconsistent. \n",
             "When alternative=\"greater\", argument \'delta\' should be positive. \n")
    }else if(alternative=="less" & delta>0){
        stop("The values given for arguments alternative and delta are inconsistent. \n",
             "When alternative=\"less\", delta should be negative. \n")
    }

    ## ** compute boundaries at decision and possibly update futility boundary at interim
    cMin <- ifelse(cNotBelowFixedc,stats::qnorm(1-alpha),-Inf)
    
    ## ## ** remove boundaries corresponding to stage that will not be reached

    if(method==1){
    
        delayedBnds <- Method1(rho_alpha = rho_alpha,
                               rho_beta = rho_beta,
                               alpha = alpha,
                               beta = beta, 
                               Kmax = kMax,
                               Info.max = NULL,
                               InfoR.i = InfoR.i,
                               InfoR.d = InfoR.d,
                               delta = delta, 
                               alternative = alternative,
                               binding=bindingFutility,
                               Trace = trace,
                               cMin = cMin)

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
                               binding=bindingFutility,
                               Trace = trace,
                               cMin = cMin)

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
                               binding=bindingFutility,
                               Trace = trace) 

    }

    ## ** output
    out <- list(call = call,
                n.obs = n,
                stage = data.frame(k = 0, type = "planning"),
                conclusion = matrix(as.character(NA), nrow = 4, ncol = kMax, dimnames = list(c("interim","reason.interim","decision","comment.decision"),NULL)),
                uk = rep(NA, kMax), 
                lk = rep(NA, kMax),
                ck = rep(NA, kMax-1),
                Info.i = rep(NA, kMax), 
                Info.d = rep(NA, kMax-1),                
                lmm = vector(mode = "list", length = kMax),
                alpha = alpha,
                alphaSpent = rep(NA, kMax),
                kMax = kMax,
                beta = beta,
                betaSpent = rep(NA, kMax),
                method = method,
                delta = data.frame(estimate = rep(NA,kMax), se = rep(NA,kMax), statistic = rep(NA,kMax), df = rep(NA,kMax), p.value = rep(NA,kMax)),
                cNotBelowFixedc = cNotBelowFixedc,
                cMin = cMin,
                bindingFutility = bindingFutility,
                alternative = alternative,
                planned = list(rho_alpha = rho_alpha,
                               rho_beta = rho_beta,
                               InflationFactor = delayedBnds$coef,
                               Info.max = delayedBnds$Info.max,
                               uk = delayedBnds$boundaries$u.k,
                               lk = delayedBnds$boundaries$l.k,
                               ck = delayedBnds$boundaries$c.k[1:(kMax-1)],
                               alphaSpent = cumsum(delayedBnds$boundaries$Inc.Type.I),
                               betaSpent = cumsum(delayedBnds$boundaries$Inc.Type.II),
                               Info.i = InfoR.i*delayedBnds$Info.max,
                               Info.d = delayedBnds$Info.d,
                               delta = delta))
    class(out) <- append("delayedGSD",class(out))
    return(out)
}


