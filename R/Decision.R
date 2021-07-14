## * Decision (documentation)
#' @title Evaluate decision at an interim, decision or final analysis of a group sequential design with delayed endpoints
#' @description Maps the statistical results at an interim, decision or final analysis into a decision regarding whether to stop recruitment (interim) or whether to reject the null hypothesis (decision/final). Stopping boundaries are updated based on observed information and correct p-values, confidence intervals and point estimates are given.
#' 
#' @param lmm object of class lmmGSD from AnalyzeData
#' @param boundaries object of class delayedGSD from CalcBoundaries
#' @param k the stage at which the decision is to be made
#' @param analysis is it an interim, decision or final analysis
#' @param Info.i the observed (where possible) or expected information at each interim and the final analysis
#' @param InfoR.d the expected or observed information ratio at each decision analysis
#' @param PositiveIsGood whether a positive effect is considered beneficial (TRUE/FALSE)
#' @param Trace whether to print some messages
#' @param bindingFutility whether to use binding futility boundaries (use TRUE for binding)
#' @param plot whether the updated boundaries and the result should be plotted
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
#' set.seed(10)
#' theData <- GenData(n=theN*2,delta=theDelta*0.8,ar=5)  #generate data with all data for in case trial completes
#' 
#' theAR <- 10  #accrual rate (pt per month)
#' theDelay <- 0.7500001  #time in months to process data
#' tau.i <- theData$d$t3[theN + ceiling(theAR*theDelay)] #time point at which to do IA
#'
#' theObsData <- SelectData(theData$d, t = tau.i, Delta.t = theDelay)  #data at IA when deciding whether to continue recruitment
#'
#' #### Analyse data at interim ####
#' myBound1 <- update(myBound0, data = theObsData, k = 1, analysis = "interim")
#' myInterim1 <- Decision(myBound1) 


## * Decision (documentation)
Decision <- function(object,
                     PositiveIsGood=TRUE, # whether positive effect is good (i.e. positive trial)
                     trace=TRUE){
  
    ## ** identify stage of the trial
    kMax <- object$kMax
    k <- object$stage$k
    if(k==kMax){
        stage <- "final"
    }else if(object$stage$decision>0){
        stage <- "decision"
    }else{
        stage  <- "interim"
    }
    
    if(trace){
        if(stage == "final"){
            cat("Decision to be take at the final analysis  (stage ",k,"). \n",sep="")
        }else if(stage == "decision"){
            cat("Decision to be take at the interim analysis of stage ",k,". \n",sep="")
        }else if(stage == "interim"){
            cat("Decision to be take at the decision analysis of stage ",k,". \n",sep="")
        }
    }

    
    ## ** extract information, test statistic, and boundaries
    ls.info <- getInformation(object, planned = FALSE)
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
    if(stage=="interim"){
        
        if(k>1 && (Info.i[k] < Info.i[k-1])){ ## skip interim analysis if information is decreasing
            object$conclusion["interim",k] <- "continue"
            object$conclusion["reason.interim",k] <- "decreasing information"
            
        } else if(Info.i[k] > Info.max){ ## skip interim analysis and continue straight to decision analysis if it is anticipated that Imax will be reached at decision analysis k
            object$conclusion["interim",k] <- "stop"
            object$conclusion["reason.interim",k] <- "Imax reached"
        } else {
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
    if(stage=="decision"){
                                        
        if(Info.d[k] < Info.i[k]){ ## if information has decreased since interim analysis
            stop("cannot handle case Info.d[k] < Info.i[k]")
        }else if(object$conclusion["reason.interim",k]=="Imax reached"){ ## if Imax has been reached at the interim
            if(Z > ck[k]){
                object$conclusion["decision",k] <- "Efficacy"
            } else {
                object$conclusion["decision",k] <- "Futility"
            }
        }else if(Info.d[k] >= Info.max){ ##  if Imax has been reached between interim and decision
            stop("don't know what to do if Id > Imax but Ii < Imax")
        } else { ## usual case
            if(Z > ck[k]){
                object$conclusion["decision",k] <- "Efficacy"
            } else {
                object$conclusion["decision",k] <- "Futility"
            }
        }
    }  

    ## ** decision at the final analysis  
    if(stage == "final"){
        if(Z > uk[k]){
            object$conclusion["decision",k] <- "Efficacy"
        } else {
            object$conclusion["decision",k] <- "Futility"
        }
    }

    ## ** export  
    return(object)
}

