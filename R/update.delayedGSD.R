## * update.delayedGSD (documentation)
#' @title Perform Next Analysis Stage in a GSD
#' @description Use the newly observed data to perform the next stage of the analysis:
#' fit a mixed model, update the information, recompute the boundaries, take a decision, and possibly correct the estimated treatment effect.
#' 
#' @param object object of class \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param data [data.frame] dataset used to fit the linear mixed model.
#' @param PositiveIsGood [logical] whether a positive effect is considered beneficial.
#' @param ddf [character] method used to compute the degrees of freedom of the test statistics. Argument passed to \code{\link{analyzeData}}.
#' @param k [integer] index of the analysis.
#' @param type.k [character] type of analysis: \code{"interim"} (after continuing recruitment),
#' \code{"decision"} (after stopping recruitment for efficacy or futility),
#' or \code{"final"} (after reaching the last stage of the trial).
#' @param trace [logical] should the execution of the function be traced?

## * update.delayedGSD (examples)
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
#' myInterim1 <- update(myBound0, data = theInterimData, k = 1, analysis = "interim")
#' print(myInterim1)
#' print(myInterim1, planned = FALSE)
#' print(myInterim1, planned = "only")
#'
#' par(mfrow = c(1,2))
#' plot(myInterim1, planned = "only")
#' plot(myInterim1)
#' 
#' #### Analyse data at the final stage ####
#' theFinalData <- SelectData(theData$d, t = 1e7, Delta.t = theDelay) 
#' myFinal <- update(myInterim1, data = theFinalData, k = 2, analysis = "final")
#' myFinal
#' print(myFinal, abreviated = FALSE)

## * update.delayedGSD (code)
#' @export
update.delayedGSD <- function(object, data, PositiveIsGood=NULL, ddf = NULL, k = NULL, type.k = NULL, trace = TRUE, ...){

    kMax <- object$kMax
    object.type.k <- object$stage$type
    object.k <- object$stage$k
    
    if(trace>0){
        if(object.type.k=="planning"){
            cat("Group Sequential Trial at the ",object.type.k," stage. \n", sep = "")
        }else{
            cat("Group Sequential Trial at the ",object.type.k," analysis of stage ",object.k,". \n", sep = "")
        }
    }

    ## ** normalize user input
    resStage <- .getStage(object.stage = object$stage,
                          object.conclusion = object$conclusion,
                          kMax = kMax,
                          k = k,
                          type.k = type.k)
    k <- resStage$k
    type.k <- resStage$type
    if(trace>0){cat("Update for the ",type.k," analysis of stage ",k,". \n", sep = "")}

    ## ** update mixed model
    if(trace>0){cat(" - fit mixed model: ", sep = "")}
    if(is.null(ddf)){
        lmm <- analyzeData(data, getinfo = TRUE, trace = trace-1)
    }else{
        lmm <- analyzeData(data, ddf = ddf, getinfo = TRUE, trace = trace-1)
    }
    if(type.k %in% c("interim","final")){
        object$lmm[[k]] <- lmm
    }else{
        object$lmm[[k+1]] <- lmm
    }
    if(trace>0){cat("done \n", sep = "")}

    ## ** update information
    if(trace>0){cat(" - udpate information: ", sep = "")}
    object <- updateInformation(object, lmm = lmm, k = k, type.k = type.k, update.stage = FALSE)
    if(trace>0){cat("done \n", sep = "")}
    
    ## ** update boundaries
    if(trace>0){cat(" - udpate boundaries: ", sep = "")}
    object <- updateBoundaries(object, k = k, type.k = type.k, trace = trace-1, update.stage = TRUE)
    if(trace>0){cat("done \n", sep = "")}

    ## ** decision
    if(trace>0){cat(" - udpate decision: ", sep = "")}
    if(!is.null(PositiveIsGood)){
        object <- Decision(object, PositiveIsGood = PositiveIsGood, k = k, type.k = type.k, trace = trace-1)
    }else{
        object <- Decision(object, k = k, type.k = type.k, trace = trace-1)
    }
    if(trace>0){cat("done \n", sep = "")}

    ## ** estimate
    if(type.k %in% c("decision","final")){
        if(trace>0){cat(" - correct estimate: ", sep = "")}

        ## FinalPvalue
        ## FinalCI
        ## FinalEstimate

        if(trace>0){cat("done \n", sep = "")}
    }
    
    ## ** export
    return(object)
}

## * .getStage
.getStage <- function(object.stage, object.conclusion, kMax, k, type.k){

    ## ** extract information
    object.k <- object.stage$k
    object.type.k <- object.stage$type
    object.conclusionInterim <- object.conclusion["interim",object.k]

    ## ** basic checks
    if(!is.null(type.k)){
        type.k <- match.arg(type.k, c("interim","decision","final"))
    }
    if(!is.null(k) && k<0 || k>kMax || k %% 1 != 0){
        stop("Argument \'k\' must be an integer between 1 and ",kMax,"")
    }

    ## ** more precise check
    if(object.k==0){
        if(!is.null(k)){
            if(k!=1){
                stop("Argument \'k\' should be 1 just after the planning stage. \n")
            }
        }else{
            k <- 1
        }
        if(!is.null(type.k)){
            if(type.k!="interim"){
                stop("Argument \'type.k\' should be \"interim\" just after the planning stage. \n")
            }
        }else{
            type.k <- "interim"
        }
    }else{
        if(is.na(object.conclusionInterim)){
            stop("Decision type.k must be made for stage ",object.k," before updating the boundaries. \n")
        }
        if(object.type.k %in% "decision"){
            stop("No more boundary to update when a decision analysis has been performed. \n")
        }
        if(object.type.k %in% "final"){
            stop("No more boundary to update when the final analysis has been performed. \n")
        }
        if(object.conclusionInterim=="continue"){
            if(!is.null(k)){
                if(k!=(object.k+1)){
                    stop("Argument \'k\' should be ",object.k+1," after continuing recruitment following the interim of stage ",object.k,". \n")
                }
            }else{
                k <- object.k + 1
            }
            if(k==kMax){
                if(!is.null(type.k)){
                    if(type.k!="final"){
                        stop("Argument \'type.k\' should be \"final\" after continuing recruitment following the interim of stage ",object.k,". \n")
                    }
                }else{
                    type.k <- "final"
                }
            }else{ ## k < kMax
                if(!is.null(type.k)){
                    if(type.k!="interim"){
                        stop("Argument \'type.k\' should be \"interim\" after continuing recruitment following the interim of stage ",object.k,". \n")
                    }
                }else{
                    type.k <- "interim"
                }
            }
        }else{ ##    object.conclusionInterim=="stop"
            if(!is.null(k)){
                if(k!=object.k){
                    stop("Argument \'k\' should be ",object.k," after stopping recruitment following the interim of stage ",object.k,". \n")
                }
            }else{
                k <- object.k
            }
            if(!is.null(type.k)){
                if(type.k!="decision"){
                    stop("Argument \'type.k\' should be \"decision\" after stopping recruitment following the interim of stage ",object.k,". \n")
                }
            }else{
                type.k <- "decision"
            }
        }
    }
    
    ## ** export
    return(list(k=k,
                type.k = type.k))
}
