#' @title Update Information in a GSD
#' @description Update the information based according to a linear mixed model.
#'
#' @param x Object of type \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param lmm Linear mixed model as a \code{lmmGSD} object, typically output from \code{\link{analyzeData}}.
#' @param k [integer] Index of the analysis.
#' @param type.k [character] Type of analysis: \code{"interim"} (after continuing recruitment),
#' \code{"decision"} (after stopping recruitment for efficacy or futility),
#' or \code{"final"} (after reaching the last stage of the trial).
#' @param update.stage [logical] should the arguments \code{k} and \code{type.k} be used to update to stage of the trial?
#' 
updateInformation <- function(object, lmm, k, type.k, update.stage = TRUE){

    ## ** update information
    if(type.k %in% c("interim","final")){
        object$Info.i[k] <- as.double(lmm$getInformation["interim"])
    }
    if(type.k %in% c("interim","decision")){
        ## update decision (even when doing interim) to ensure monotone information
        object$Info.d[k] <- as.double(lmm$getInformation["decision"])
    }

    ## ** udpate.stage
    if(update.stage){
        object$stage$k <- k
        object$stage$type <- type.k
    }

    ## ** export
    return(object)
}
