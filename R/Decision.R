## * Decision (documentation)
#' @title Take Decision in a GSD
#' @description Maps the statistical results at an interim, decision or final analysis into a decision regarding whether to stop recruitment (interim) or whether to reject the null hypothesis (decision/final). Stopping boundaries are updated based on observed information and correct p-values, confidence intervals and point estimates are given.
#' 
#' @param object object of class \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param PositiveIsGood [logical] whether a positive effect is considered beneficial (TRUE/FALSE)
#' @param k [integer] index of the analysis.
#' @param type.k [character] type of analysis: \code{"interim"} (after continuing recruitment),
#' \code{"decision"} (after stopping recruitment for efficacy or futility),
#' or \code{"final"} (after reaching the last stage of the trial).
#' @param trace [logical] should the execution of the function be traced?
#' 
#' @details Function that maps the statistical results at interim (or final) analysis into decision to reject, or continue or stop inclusind subjects.
#' @return ff
#' @author Paul Blanche
#' 
## * Decision (example)
#' @examples
#'
#' #### Planning #####
#' theAlpha <- 0.025
#' theBeta <- 0.2
#' theDelta <- 1.5
#' theK <- 2
#' theN <- 82
#' 
#' myBound0 <- CalcBoundaries(kMax=theK,
#'                      sided=1,
#'                      alpha=theAlpha,
#'                      beta=theBeta,
#'                      InfoR.i=c(0.5,1),
#'                      gammaA=2,
#'                      gammaB=2,
#'                      method=1,
#'                      delta=theDelta,
#'                      InfoR.d=0.55)
#'
#' #### Simulate data ####
#' ## generate data with all data for in case trial completes
#' set.seed(10)
#' theData <- GenData(n=theN*2,delta=theDelta*0.8,ar=5)  
#' 
#' theAR <- 10  #accrual rate (pt per month)
#' theDelay <- 0.7500001  #time in months to process data
#' tau.i <- theData$d$t3[theN + ceiling(theAR*theDelay)] #time point at which to do IA
#'
#' theObsData <- SelectData(theData$d, t = tau.i, Delta.t = theDelay)
#'
#' #### Analyse data at interim ####
#' myLMM <- analyzeData(theObsData)
#' myBound1 <- updateBoundaries(myBound0, lmm = myLMM, k = 1, type.k = "interim")
#' myInterim1 <- Decision(myBound1, k = 1, type.k = "interim")
#' myInterim1$conclusion

## * Decision (documentation)
#' @export
Decision <- function(object,
                     k,
                     type.k,
                     PositiveIsGood = TRUE, # whether positive effect is good (i.e. positive trial)
                     trace = TRUE){
  
    ## ** identify stage of the trial
    kMax <- object$kMax
    
    if(trace){
        if(type.k == "final"){
            cat("Decision to be taken at the final analysis  (stage ",k,") \n",sep="")
        }else if(type.k == "decision"){
            cat("Decision to be taken at the decision analysis of stage ",k," \n",sep="")
        }else if(type.k == "interim"){
            cat("Decision to be taken at the interim analysis of stage ",k," \n",sep="")
        }
    }

    
    ## ** extract information, test statistic, and boundaries
    ls.info <- getInformation(object, planned = TRUE)
    Info.max <- ls.info$Info.max
    Info.i <- ls.info$Info.i
    Info.d <- ls.info$Info.d
    uk <- ls.info$uk
    lk <- ls.info$lk
    ck <- ls.info$ck
    
    Z <- ls.info$delta[NROW(ls.info$delta),"statistic"]

    if(!PositiveIsGood){
        Z <- -Z
        if(trace){message("Negative effect is good: change sign of Z")}
    }

    ## ** decision at interim
    if(type.k=="interim"){
        
        if(k>1 && (Info.i[k] < Info.i[k-1])){ ## will continue to the next interim when information decreased
            object$conclusion["interim",k] <- "continue"
            object$conclusion["reason.interim",k] <- "decreasing information"
        } else if((Info.i[k] >= Info.max) || (Info.d[k] >= Info.max)){ ## stop and do decision when Imax is already reached or anticipated to be reached
            object$conclusion["interim",k] <- "stop"
            object$conclusion["reason.interim",k] <- "Imax reached"
        } else { ## usual case
            if(Z>uk[k]){
                object$conclusion["interim",k] <- "stop"
                object$conclusion["reason.interim",k] <- "efficacy"
            } else if (Z<lk[k]){
                object$conclusion["interim",k] <- "stop"
                object$conclusion["reason.interim",k] <- "futility"
            } else {
                object$conclusion["interim",k] <- "continue"
                object$conclusion["reason.interim",k] <- "no boundary crossed"
            }
        }
    }
    
    ## ** decision at decision analysis
    if(type.k=="decision"){
                                        
        if(Info.d[k] < Info.i[k]){ ## if information has decreased since interim analysis
            stop("cannot handle case Info.d[k] < Info.i[k]")
        }else if((object$conclusion["reason.interim",k]=="Imax reached") || (Info.d[k] >= Info.max)){ ## if Imax has been reached at the interim
            if(Z > ck[k]){
                object$conclusion["decision",k] <- "Efficacy"
            } else {
                object$conclusion["decision",k] <- "Futility"
            }
        } else { ## usual case
            if(Z > ck[k]){
                object$conclusion["decision",k] <- "Efficacy"
            } else {
                object$conclusion["decision",k] <- "Futility"
            }
        }
    }  

    ## ** decision at the final analysis  
    if(type.k == "final"){
        if(Z > uk[k]){
            object$conclusion["decision",k] <- "Efficacy"
        } else {
            object$conclusion["decision",k] <- "Futility"
        }
    }

    ## ** export  
    return(object)
}

