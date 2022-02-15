#' @title Update Information in a GSD
#' @description Update the information based according to a linear mixed model.
#'
#' @param object Object of type \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param lmm Linear mixed model as a \code{lmmGSD} object, typically output from \code{\link{analyzeData}}.
#' @param k [integer] Index of the analysis.
#' @param n.decision [integer] expected number of patients at decision analysis.
#' @param type.k [character] Type of analysis: \code{"interim"} (after continuing recruitment),
#' \code{"decision"} (after stopping recruitment for efficacy or futility),
#' or \code{"final"} (after reaching the last stage of the trial).
#' @param update.stage [logical] should the arguments \code{k} and \code{type.k} be used to update to stage of the trial?
#' 
updateInformation <- function(object, lmm, n.decision, k, type.k, update.stage = TRUE){

    kMax <- object$kMax
    
    ## ** update information at current stage
    if(type.k %in% c("interim","final")){
        object$Info.i[k] <- as.double(lmm$getInformation["interim"])
    }
    if(type.k %in% c("interim","decision")){
        ## update decision (even when doing interim) to ensure monotone information
        if(is.null(n.decision)){
            object$Info.d[k] <- as.double(lmm$getInformation["decision"])
        }else{
            
            object$Info.d[k] <- (n.decision/lmm$n["total"])*as.double(lmm$getInformation["decision"])
        }
    }

    ## ## ** update information at future stages
    ## Has some side effect: change the boundary at the current stage if information exceed Imax
    ## if(type.k == "interim"){
        
    ##     if((k+1)<=(kMax-1)){
    ##         object$Info.i[(k+1):(kMax-1)] <- object$Info.i / (object$planned$Info.i[k]/object$planned$Info.i[(k+1):(kMax-1)])
    ##         object$Info.d[(k+1):(kMax-1)] <- object$Info.d / (object$planned$Info.d[k]/object$planned$Info.d[(k+1):(kMax-1)])
    ##     }
    ##     object$Info.i[kMax] <- object$Info.d / (object$planned$Info.d[k]/object$planned$Info.i[kMax])
        
    ## }else if(type.k == "decision"){ ## will never reach any other stage

    ##     object$Info.i[k:kMax] <- NA
    ##     if((k+1)<=(kMax-1)){
    ##         object$Info.d[(k+1):(kMax-1)] <- NA
    ##     }

    ## }

    ## ** udpate.stage
    if(update.stage){
        object$stage$k <- k
        object$stage$type <- type.k
    }

    ## ** export
    return(object)
}
