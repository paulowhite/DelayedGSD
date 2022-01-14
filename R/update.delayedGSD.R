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
#' @param p.value [logical] should the p-value be computed at decision?
#' @param ci [logical] should the confidence intervalsbe computed at decision?
#' @param estimate [logical] should a de-biased estimate be computed at decision? WARNING: this is experiment and not reliable.
#' @param trace [logical] should the execution of the function be traced?
#' @param ... not used, for compatibility with the generic method.

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
#' 
#' myBound0 <- CalcBoundaries(kMax=theK,
#'                      alpha=theAlpha,
#'                      beta=theBeta,
#'                      InfoR.i=c(0.5,1),
#'                      rho_alpha=2,
#'                      rho_beta=2,
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
#' myInterim1 <- update(myBound0, data = theInterimData) ## k = 1, analysis = "interim"
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
#' myFinal <- update(myInterim1, data = theFinalData) ## k = 2, analysis = "final"
#' myFinal
#' print(myFinal, abreviated = FALSE)
#' plot(myFinal)

## * update.delayedGSD (code)
#' @export
update.delayedGSD <- function(object, data, PositiveIsGood=NULL, ddf = NULL, k = NULL, type.k = NULL,
                              p.value = TRUE, ci = TRUE, estimate = TRUE,
                              trace = TRUE, ...){

    

    ## ** normalize user input
    kMax <- object$kMax
    object.type.k <- object$stage$type
    object.k <- object$stage$k

    if(object.type.k %in% "decision"){
        stop("No more boundary to update when a decision analysis has been performed. \n")
    }
    if(object.type.k %in% "final"){
        stop("No more boundary to update when the final analysis has been performed. \n")
    }
    resStage <- .getStage(object.stage = object$stage,
                          object.conclusion = object$conclusion,
                          kMax = kMax,
                          k = k,
                          type.k = type.k,
                          nextStage = TRUE)
    k <- resStage$k
    type.k <- resStage$type

    if(trace>0){
        if(object.type.k=="planning"){
            cat("Group Sequential Trial at the ",object.type.k," stage. \n", sep = "")
        }else{
            cat("Group Sequential Trial at the ",object.type.k," analysis of stage ",object.k,". \n", sep = "")
        }
        cat("Update for the ",type.k," analysis of stage ",k,". \n", sep = "")
    }

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
    if(trace>0){cat(" - update information: ", sep = "")}
    object <- updateInformation(object, lmm = lmm, k = k, type.k = type.k, update.stage = FALSE)
    if(trace>0){cat("done \n", sep = "")}

    ## ** update boundaries
    if(trace>0){cat(" - update boundaries: ", sep = "")}
    object <- updateBoundaries(object, k = k, type.k = type.k, trace = trace-1, update.stage = TRUE)
    if(trace>0){cat("done \n", sep = "")}

    ## ** decision
    if(trace>0){cat(" - update decision: ", sep = "")}
    if(type.k != "interim" || is.na(object$conclusion["interim",k])){
        if(!is.null(PositiveIsGood)){
            object <- Decision(object, PositiveIsGood = PositiveIsGood, k = k, type.k = type.k, trace = trace-1)
        }else{
            object <- Decision(object, k = k, type.k = type.k, trace = trace-1)
        }
    }
    
    if(trace>0){cat("done \n", sep = "")}

    ## ** estimate
    if(type.k %in% c("decision","final")){
        if(trace>0){cat(" - correct estimate: ", sep = "")}

        delta <- getInformation(object)$delta

        object$correction <- data.frame(estimate=NA,
                                        lower=NA,
                                        upper=NA,
                                        p.value=NA)

        ## *** p.value
        if(p.value){
            object$correction$p.value <- FinalPvalue(Info.d = object$Info.d,  
                                                     Info.i = object$Info.i,  
                                                     ck = object$ck,   
                                                     lk = object$kk,  
                                                     uk = object$uk,  
                                                     kMax = kMax, 
                                                     delta = 0,  
                                                     estimate = delta[NROW(delta),"estimate"])
        }
        
        ## *** CI
        if(ci){
            resCI <- FinalCI(Info.d = object$Info.d,  
                             Info.i = object$Info.i,  
                             ck = object$ck,   
                             lk = object$kk,  
                             uk = object$uk,  
                             kMax = kMax, 
                             alpha = object$alpha,  
                             estimate = delta[NROW(delta),"estimate"])
            object$correction$lower <- resCI["lower"]
            attr(object$correction$lower,"error") <- attr(resCI,"error")["lower"]
            object$correction$upper <- resCI["upper"]
            attr(object$correction$lower,"error") <- attr(resCI,"error")["upper"]
        }
        
        ## *** Estimate
        if(estimate){
            object$correction$estimate <- FinalEstimate(Info.d = object$Info.d,  
                                                        Info.i = object$Info.i,  
                                                        ck = object$ck,   
                                                        lk = object$kk,  
                                                        uk = object$uk,  
                                                        kMax = kMax, 
                                                        estimate = delta[NROW(delta),"estimate"])
        }
        
        if(trace>0){cat("done \n", sep = "")}
    }
    
    ## ** export
    return(object)
}

## * .getStage
## Identify the current stage or check that the values of k and type.k are compatible with the object
.getStage <- function(object.stage, object.conclusion, kMax, k, type.k, nextStage){

    ## ** extract information
    object.k <- object.stage$k
    object.type.k <- object.stage$type

    ## ** basic checks
    if(!is.null(type.k)){
        type.k <- match.arg(type.k, c("interim","decision","final"))
    }
    if(!is.null(k)){
        if(k<0 || k>kMax || k %% 1 != 0){
            stop("Argument \'k\' must be an integer between 1 and ",kMax,"")
        }
    }

    ## ** more precise check
    if(object.k==0){

        ## *** planning
        if(nextStage){
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
            if(!is.null(k)){
                if(k!=0){
                    stop("Argument \'k\' should be 0 at the planning stage. \n")
                }
            }else{ 
                k <- 0
            }
            if(!is.null(type.k)){
                if(type.k!="planning"){
                    stop("Argument \'type.k\' should be \"planning\" just at the planning stage. \n")
                }
            }else{
                type.k <- "planning"
            }
        }
    }else if(object.type.k=="interim"){

        if(nextStage){
            object.conclusionInterim <- object.conclusion["interim",object.k]

            if(object.conclusionInterim=="continue"){ ## at interim where we conclude to continue
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
            }else if(object.conclusionInterim=="stop"){ ## at interim where we conclude to stop
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

        }else{
            if(!is.null(k)){
                if(k %in% 1:object.k == FALSE){
                    stop("Argument \'k\' should be an integer between 1 and ",object.k,". \n")
                }
            }else{
                k <- object.k
            }
            if(!is.null(type.k)){
                if(type.k %in% "interim" == FALSE){
                    stop("Argument \'type.k\' can only be \"interim\". \n")
                }
            }else{
                type.k <- "interim"
            }
        }

    }else if(object.type.k=="decision"){
        if(nextStage){
            stop("The decision analysis has already been reached, no further analysis. \n")
        }
        if(!is.null(k)){
            if(k %in% 1:object.k == FALSE){
                stop("Argument \'k\' should be an integer between 1 and ",object.k,". \n")
            }
        }else{
            k <- object.k
        }
        if(!is.null(type.k)){
            if(k<object.k && type.k=="decision"){
                stop("Argument \'type.k\' cannot be \"decision\" as recruitement continued after interim ",k,". \n")
            }else if(type.k %in% c("decision","interim") == FALSE){
                stop("Argument \'type.k\' can either be \"interim\" or \"decision\". \n")
            }
        }else{
            type.k <- "decision"
        }

    }else if(object.type.k=="final"){
    
        if(nextStage){
            stop("The final analysis has already been reachedm no further analysis. \n")
        }
        if(!is.null(k)){
            if(k %in% 1:object.k == FALSE){
                stop("Argument \'k\' should be an integer between 1 and ",object.k,". \n")
            }
        }else{
            k <- object.k
        }
        if(!is.null(type.k)){
            if(k==object.k){
                if(type.k!="final"){
                    stop("Argument \'type.k\' must be \"final\" at stage ",k,". \n")
                }
            }else if(type.k %in% c("decision","interim") == FALSE){
                if(k>2){
                    stop("Argument \'type.k\' can either be \"interim\" or \"decision\" between stage 1 and ",k-1,". \n")
                }else{
                    stop("Argument \'type.k\' can either be \"interim\" or \"decision\" at stage 1. \n")
                }
            }
        }else{
            type.k <- "final"
        }

    }
    
    ## ** export
    return(list(k=k,
                type.k = type.k))
}
