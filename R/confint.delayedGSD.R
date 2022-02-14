## * confint.delayedGSD (documentation)
#' @title Treatment Effect Estimates for a Group Sequential Design with Delayed Endpoints
#' @description Extract estimate, confidence intervals, p-value, ... relative to the treatment effect at a specific stage of a group sequential design with delayed endpoints.
#' By default extract value for the latest stage that has been performed.
#' 
#' @param object object of class \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param parm not used, for compatibility with the generic method.
#' @param level not used, for compatibility with the generic method.
#' @param method [character] type of estimate to output: can be \code{"ML"} or  \code{"corrected ML"}, the latter accounting for the group sequential design.
#' @param k [integer] stage relative to which the estimates should be output.
#' @param type.k [character] type  of stage relative to which the estimates should be output: \code{"interim"}, \code{"decision"}, or \code{"final"}.
#' @param ... not used, for compatibility with the generic method.
#' 

## * confint.delayedGSD (code)
#' @export
confint.delayedGSD <- function(object, parm = NULL, level = NULL, method = NULL, k = NULL, type.k = NULL, ...){


    ## ** check user input
    if(object$stage$type=="planning"){return(NULL)} ## no estimate at the planning stage
    if(!is.null(parm)){
        stop("Argument \'parm\' not used. ")
    }
    if(!is.null(level)){
        stop("Argument \'level\' not used. ")
    }
    dots <- list(...)
    if(length(dots)){
        stop("Argument(s) \'",paste(names(dots), collapse = "\' \'"),"\' not used. ")
    }

    kMax <- object$kMax
    resStage <- .getStage(object.stage = object$stage,
                          object.conclusion = object$conclusion,
                          kMax = kMax,
                          k = k,
                          type.k = type.k,
                          nextStage = FALSE)

    if(is.null(k)){
        k <- resStage$k
    }
    if(is.null(type.k)){
        type.k <- resStage$type.k
    }

    if(!is.null(method)){
        method <- match.arg(method, c("ML","corrected ML"), several.ok = TRUE)
        ## if("corrected ML" %in% method &&  type.k == "interim"){
        ##     stop("Argument \'method\' cannot be \"corrected ML\" at interim. \n")
        ## }
    }

    ## ** extract information from object
    if(type.k == "decision"){
        iDelta <- object$delta[k+1,]
    }else{
        iDelta <- object$delta[k,]
    }
    out <- data.frame(method = "ML",
                      stage = k,
                      type = type.k,
                      coef = object$lmm[[1]]$name.coef,
                      estimate = iDelta$estimate,
                      se = iDelta$se,
                      lower = iDelta$estimate + stats::qt(object$alpha/2, df = iDelta$df) * iDelta$se,
                      upper = iDelta$estimate + stats::qt(1-object$alpha/2, df = iDelta$df) * iDelta$se,
                      statistic = iDelta$statistic,
                      df = iDelta$df,
                      p.value = iDelta$p.value)
    if(!is.null(object$correction)){
        out <- rbind(out, data.frame(method = "corrected ML",
                                     stage = k,
                                     type = type.k,
                                     coef = object$lmm[[1]]$name.coef,
                                     estimate = object$correction$estimate,
                                     se = NA,
                                     lower = object$correction$lower,
                                     upper = object$correction$upper,
                                     statistic = iDelta$statistic,
                                     df = NA,
                                     p.value = object$correction$p.value)
                     )

    }

    ## ** export
    if(is.null(method)){
        out <- out[NROW(out),,drop=FALSE]
    }else{
        out <- out[out$method %in% method,,drop=FALSE]
    }
    rownames(out) <- NULL            
    return(out)

}
