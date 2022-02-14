## * coef.delayedGSD (documentation)
#' @title Treatment Effect Estimates for a Group Sequential Design with Delayed Endpoints
#' @description Extract estimate relative to the treatment effect at a specific stage of a group sequential design with delayed endpoints.
#' By default extract value for the latest stage that has been performed.
#' 
#' @param object object of class \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param type [character] Shoudl the estimate effect (\code{"effect"}), boundaries (\code{"boundaries"}), information  (\code{"information"}), or decision  (\code{"deicsion"})
#' be output. The estimate is only displayed for the current stage while the other are displayed for all stages.
#' @param planned [character] Should only the boundary/information used to plan the trial be output (\code{"only"}),
#' or only the estimated information (\code{FALSE}),
#' or the estimated information when available and otherwise the information used to plan the trial.
#' @param k [integer] stage relative to which the estimates should be output.
#' @param type.k [character] type of stage relative to which the estimates should be output: \code{"interim"}, \code{"decision"}, or \code{"final"}.
#' @param method [character] type of estimate to output: can be \code{"ML"} or  \code{"corrected ML"}, the latter accounting for the group sequential design.
#' @param ... not used, for compatibility with the generic method.
#' 

## * coef.delayedGSD
##' @export
coef.delayedGSD <- function(object, type = "effect", planned = NULL, k = NULL, type.k = NULL, method = NULL, ...){

    ## ** extract information from object
    type <- match.arg(type, c("effect","boundary","information","decision"))
    
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

    test.planning <- (type.k=="planning") || identical(planned,"only")
    if(is.null(planned)){
        planned <- test.planning
    }
    Info.max <- object$planned$Info.max

    ## ** retrive quantities
    if(type=="effect"){
        outConfint <- stats::confint(object, method = method, k = k, type.k = type.k)
        out <- stats::setNames(outConfint$estimate,outConfint$method)
    }else if(type=="boundary"){
        if(test.planning){

            out <- data.frame(Stage = 1:kMax,
                              Fbound = object$planned$lk,
                              Ebound = object$planned$uk,
                              Cbound = c(object$planned$ck,utils::tail(object$planned$uk,1)))

        }else{
            statistic <- unname(object$delta$statistic)
            statistic.interim <- c(statistic[1:k], rep(NA,kMax-k))
            statistic.interim[kMax] <- NA ## at final analysis no interim only decision
            statistic.decision <- rep(NA,kMax)
            if(type.k=="decision"){
                statistic.decision[k] <- statistic[k+1]
            }else if(type.k=="final"){
                statistic.decision[k] <- statistic[k]
            }

            ## remove critical boundaries when continuing at interim
            if(planned==FALSE){
                object$ck[which(object$conclusion["interim",]=="continue")] <- NA
            }        
                
            ## assemble
            out <- data.frame(Stage = 1:kMax, Fbound = object$lk, Ebound = object$uk, statistic.interim = statistic.interim,
                              Cbound = c(object$ck, utils::tail(object$uk,1)), statistic.decision = statistic.decision)
            rownames(out) <- NULL                
                
            }

    }else if(type=="information"){

        if(test.planning){
            Info.i <- object$planned$Info.i ## interim and final
            Info.d <- c(object$planned$Info.d,utils::tail(object$planned$Info.i,1)) ## decision and final
            InfoR.i <- Info.i/Info.max
            InfoR.d <- Info.d/Info.max
            if(!is.null(object$n.obs)){
                n.obs <- object$n.obs*InfoR.d
            }else{
                n.obs <- NULL
            }
        }else{
            Info.i <- object$Info.i ## interim and final
            Info.d <- c(object$Info.d, utils::tail(object$Info.i,1)) ## decision and final
            InfoR.i <- Info.i/Info.max
            InfoR.d <- Info.d/Info.max
            n.obs <- unlist(lapply(object$lmm, function(iLMM){
                if(!is.null(iLMM)){as.double(iLMM$sample.size["decision"])}else{NA}
            }))
            index.lastNNA <- utils::tail(which(!is.na(n.obs)),1)
            index.NA <- which(is.na(n.obs))
            n.obs[index.NA] <- n.obs[index.lastNNA] * InfoR.d[index.NA]/InfoR.d[index.lastNNA]
        }
        out <- data.frame(Stage = 1:kMax,
                          Interim = Info.i,
                          Interim.pc = InfoR.i,
                          Decision = Info.d,
                          Decision.pc = InfoR.d)
        if(!is.null(n.obs)){
            out  <- cbind(out, n = n.obs)
        }
        rownames(out) <- NULL
        attr(out,"Info.max") <- Info.max
    }else if(type == "decision"){
        out <- stats::setNames(object$conclusion["interim",], paste0("stage ",1:kMax))
        if(any(!is.na(object$conclusion["decision",]))){
            toadd <- object$conclusion["decision",]
            out[!is.na(toadd)] <- toadd[!is.na(toadd)]
        }
    }

    ## ** export
    return(out)
}
