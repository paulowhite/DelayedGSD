## * updateBoundaries (documentation)
#' @title Update Boundaries of a GSD
#' @description Recompute the boundaries based according to the current information.
#'
#' @param object Object of type \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param lmm Optional argument used to update the information. Object of type \code{lmmGSD}, typically output from \code{\link{analyzeData}}.
#' @param k [integer] Index of the analysis.
#' @param type.k [character] Type of analysis: \code{"interim"} (after continuing recruitment),
#' \code{"decision"} (after stopping recruitment for efficacy or futility),
#' or \code{"final"} (after reaching the last stage of the trial).
#' @param trace [logical] should the execution of the function be traced?
#' @param update.stage [logical] should the arguments \code{k} and \code{type.k} be used to update to stage of the trial?
#' 

## * updateBoundaries (example)
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
#'                           sided=1,
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
#'
#' #### Analyse data at the first interim ####
#' theInterimData <- SelectData(theData$d, t = tau.i, Delta.t = theDelay)
#' 
#' myLMM <- analyzeData(theInterimData)
#' myBound1 <- updateBoundaries(myBound0, lmm = myLMM, k = 1, type.k = "interim")
#' print(myBound1)
#' print(myBound1, planned = FALSE)
#' print(myBound1, planned = "only")
#'
#' par(mfrow = c(1,2))
#' plot(myBound1, planned = "only")
#' plot(myBound1)
#' 
#' #### Analyse data at the final ####
#' theFinalData <- SelectData(theData$d, t = 1e7, Delta.t = theDelay) 
#' 
#' myLMM <- analyzeData(theFinalData)
#' myBound2 <- updateBoundaries(myBound1, lmm = myLMM, k = 2, type.k = "final")
#' plot(myBound2)

## * updateBoundaries (code)
#' @export
updateBoundaries <- function(object, lmm = NULL, k, type.k, update.stage = TRUE, trace = TRUE){

    kMax <- object$kMax
    Info.max <- object$Info.max
    uk <- object$uk
    lk <- object$lk
    ck <- object$ck
    bindingFutility <- object$bindingFutility
    
    ## ** [optional] update information
    if(!is.null(lmm)){
        if(type.k %in% c("interim","final")){
            object$lmm[[k]] <- lmm
        }else{
            object$lmm[[k+1]] <- lmm
        }
        ## object$Info.i [1] 2.193777 3.640132
        object <- updateInformation(object, lmm = lmm, k = k, type.k = type.k, update.stage = FALSE)
    }

    Info.i <- object$Info.i
    Info.d <- object$Info.d

    ## ** update boundaries
    ## *** at interim
    if(type.k == "interim"){

        if(k>1 && (Info.i[k] < Info.i[k-1])){ ## continue to the next interim when information decreased
            object$lk[k]  <- NA
            object$uk[k]  <- NA
            object$ck[k]  <- NA
        }else if((Info.i[k] >= Info.max) || (Info.d[k] >= Info.max)){ ## stop and do decision when Imax is already reached or anticipated to be reached
            object$lk[k]  <- NA
            object$uk[k]  <- NA
            object$ck[k]  <- NA
        }else{ ## usual evaluation

            ## ## Handling information greater than Info.max at future interim or decision analyses or at the final analysis
            ## if(any(Info.i>object$Info.max)){
            ##     attr(Info.i,"InflationFactor") <- object$InflationFactor
            ## }

            newBounds <- CalcBoundaries(kMax=kMax,
                                        sided=object$sided,
                                        alpha=object$alpha, 
                                        beta=object$beta,  
                                        InfoR.i=Info.i/object$Info.max,  
                                        gammaA=object$gammaA,
                                        gammaB=object$gammaB,
                                        method=object$method,
                                        cNotBelowFixedc=object$cNotBelowFixedc, 
                                        delta=object$delta,  
                                        InfoR.d=Info.d/object$Info.max,  
                                        trace=trace,
                                        bindingFutility = bindingFutility)

            object$lk  <- newBounds$lk
            object$uk  <- newBounds$uk
            object$ck  <- newBounds$ck
        } 

    }

    ## *** at decision
    if(type.k == "decision"){

        ## Handling information greater than Info.max at future interim or decision analyses or at the final analysis:
        ## we will never reach these analyses so no need to compute boundaries there and we should not need the information at those stages
        Info.i[(k+1):kMax] <- NA ## future interim or final
        if((k+1)>=(kMax-1)){Info.d[(k+1):(kMax-1)] <- NA} ## future decision


        
        if(Info.d[k] < Info.i[k]){ ## if information has decreased since interim analysis
            ## the real problematic case is when Info.d[k] < Imax but Info.i[k] > Imax
            ## What to do?
            stop("Do not know how to deal when Information decreases between interim and decision. \n")
        }else if((object$conclusion["reason.interim",k]=="Imax reached") || (Info.d[k] >= Info.max)){ ## if Imax has been reached at the interim and continue to increase, or has been reached at decision
            ## recompute the decision boundary to spend all the alpha
            if(object$method==1){
                newBounds2 <- Method1(uk=uk[1:k],
                                      lk=lk[1:k],
                                      Info.i=Info.i[1:k],
                                      Info.d=Info.d[k],
                                      Info.max=Info.max,
                                      sided=object$sided,
                                      ImaxAnticipated=TRUE,
                                      rho=object$gammaA,
                                      alpha=object$alpha,
                                      bindingFutility=bindingFutility)
                object$ck[k:(kMax-1)]  <- c(newBounds2[k], rep(NA,kMax-k-1))
            }else{
                stop("Method ",object$method," not implemented when Imax reached. \n")
            }
        }else { ## usual evaluation
            newBounds <- CalcBoundaries(kMax=kMax,
                                        sided=object$sided,
                                        alpha=object$alpha, 
                                        beta=object$beta,  
                                        InfoR.i=Info.i/object$Info.max,  
                                        gammaA=object$gammaA,
                                        gammaB=object$gammaB,
                                        method=object$method,
                                        cNotBelowFixedc=object$cNotBelowFixedc, 
                                        delta=object$delta,  
                                        InfoR.d=Info.d/object$Info.max,  
                                        trace=trace,
                                        bindingFutility = bindingFutility)

            object$lk  <- newBounds$lk
            object$uk  <- newBounds$uk
            object$ck  <- newBounds$ck
        }
    }
    
    ## *** at final
    if(type.k == "final"){ ##  that could be replace by a call to CalcBoundaries
        if(bindingFutility){test.type <- 3}else{test.type <- 4}

        StandardDesign <- gsDesign::gsDesign(k = k, test.type = test.type, alpha = object$alpha, beta = object$beta,
                                             timing = c(Info.i[1:(k-1)], object$Info.max)/object$Info.max,
                                             n.fix = 1, n.I = Info.i / (object$Info.max/object$InflationFactor), maxn.IPlan = object$InflationFactor,
                                             ## n.fix=object$Info.max/object$InflationFactor, n.I=Info.i, maxn.IPlan=object$Info.max ## should be equivalent
                                             sfu = gsDesign::sfPower, sfupar = object$gammaA,
                                             sfl = gsDesign::sfPower, sflpar = object$gammaB)

        object$lk[k]  <- StandardDesign$upper$bound[k]
        object$uk[k]  <- StandardDesign$upper$bound[k]
    }


    ## ** export stage
    if(update.stage){
        object$stage$k <- k
        object$stage$type <- type.k
    }

    ## ** export object
    return(object)
}
