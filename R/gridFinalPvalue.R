### plotFinalPvalue.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 21 2023 (09:42) 
## Version: 
## Last-Updated: okt 30 2023 (11:42) 
##           By: Brice Ozenne
##     Update #: 94
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * gridFinalPvalue (documentation)
##' @title P-value over the Ordering Space
##' 
##' @description Evaluate the p-value over stages and test statistics, following armitage ordering (Armitage 1957).
##' 
##' @param object 
##' @inheritParams FinalPvalue
##' 
##' @return a \code{data.frame} containing \itemize{
##' \item \code{k} the stage
##' \item \code{z} the test statistic value
##' \item \code{p.value} the p-value
##' }
##'
##' @references P. Armitage, Restricted sequential procedures. Biometrika 1957 (9-56)
##' 
##' @examples
##' ## 2 stages, no fixC
##' myBound2 <- CalcBoundaries(kMax = 2,  
##'                          alpha = 0.025, 
##'                          beta = 0.2,  
##'                          InfoR.i = c(0.55, 1.00),  
##'                          InfoR.d = c(0.62, 1),  
##'                          rho_alpha = 2,  
##'                          rho_beta = 2,  
##'                          method = 1,  
##'                          cNotBelowFixedc = FALSE,
##'                          bindingFutility= TRUE,
##'                          delta = 0.6)
##' grid.p2 <- gridFinalPvalue(myBound2)
##' grid.p2
##' 
##' ## 2 stages, fixC
##' myBound2.fixC <- CalcBoundaries(kMax = 2,  
##'                          alpha = 0.025, 
##'                          beta = 0.2,  
##'                          InfoR.i = c(0.55, 1.00),  
##'                          InfoR.d = c(0.62, 1),  
##'                          rho_alpha = 2,  
##'                          rho_beta = 2,  
##'                          method = 1,  
##'                          cNotBelowFixedc = TRUE,
##'                          bindingFutility= TRUE,
##'                          delta = 0.6)
##' gridFinalPvalue(myBound2.fixC,
##'                 continuity.correction = 0)
##' gridFinalPvalue(myBound2.fixC,
##'                 continuity.correction = 1)
##' gridFinalPvalue(myBound2.fixC,
##'                 continuity.correction = 2)
##'
##' ## 3 stages, no fixC
##' myBound3 <- CalcBoundaries(kMax = 3,  
##'                          alpha = 0.025, 
##'                          beta = 0.2,  
##'                          InfoR.i = c(0.45, 0.65, 1.00),  
##'                          InfoR.d = c(0.52, 0.72, 1),  
##'                          rho_alpha = 2,  
##'                          rho_beta = 2,  
##'                          method = 1,  
##'                          cNotBelowFixedc = FALSE,
##'                          bindingFutility= TRUE,
##'                          delta = 0.6)
##' grid.p3 <- gridFinalPvalue(myBound3)
##' grid.p3
##' 
##' ## 3 stages, fixC
##' myBound3.fixC <- CalcBoundaries(kMax = 3,  
##'                          alpha = 0.025, 
##'                          beta = 0.2,  
##'                          InfoR.i = c(0.45, 0.65, 1.00),  
##'                          InfoR.d = c(0.52, 0.72, 1),  
##'                          rho_alpha = 2,  
##'                          rho_beta = 2,  
##'                          method = 1,  
##'                          cNotBelowFixedc = TRUE,
##'                          bindingFutility= TRUE,
##'                          delta = 0.6)
##' gridFinalPvalue(myBound3.fixC, digits = 3,
##'                 continuity.correction = 0)
##' gridFinalPvalue(myBound3.fixC, digits = 3,
##'                 continuity.correction = 1)
##' gridFinalPvalue(myBound3.fixC, digits = 3,
##'                 continuity.correction = 2)
##'

## * gridFinalPvalue (code)
##' @export
gridFinalPvalue <- function(object,
                            continuity.correction,
                            seq.futility = c(-3,0),
                            seq.kMax = c(-3,1,3,6),
                            seq.efficacy = c(2.5,3.25,4,5),
                            plot = TRUE,
                            col = NULL,
                            title = NULL,
                            xlim = NULL,
                            pch = 20,
                            cex = 2,
                            lwd = 2,
                            digits = 2){

    ## ** normalize arguments
    options <- DelayedGSD.options()
    FCT.p_value <- options$FCT.p_value
    seq.futility <- sort(seq.futility)
    seq.kMax <- sort(seq.kMax)
    seq.efficacy <- sort(seq.efficacy)

    #browser()
    
    if(inherits(object,"delayedGSD")){

        kMax <- object$kMax
        method <- object$method
        bindingFutility <- object$bindingFutility
        cNotBelowFixedc <- object$cNotBelowFixedc

        if(object$stage$type=="planning"){
            kCurrent <- kMax
            Info.d <- object$planned$Info.d
            Info.i <- object$planned$Info.i
            ck <- object$planned$ck
            ck.unrestricted <- object$planned$ck.unrestricted
            lk <- object$planned$lk
            uk <- object$planned$uk
            reason.interim <- rep("",kMax)
            alphaSpent <- object$planned$alphaSpent
        }else{
            kCurrent <- object$stage["k"]
            if(kCurrent!=kMax){
                stop("Current stage (k=",kCurrent,") is not kMax=",kMax,". \n",
                     "Can only display p-values over the whole space. \n")
            }
            Info.d <- object$Info.d
            Info.i <- object$Info.i
            ck <- object$ck
            ck.unrestricted <- object$ck.unrestricted
            lk <- object$lk
            uk <- object$uk
            reason.interim <- object$conclusion["reason.interim",]
            alphaSpent <- object$alphaSpent
        }
    }else if(is.list(object)){

        valid.args <- c("Info.d", "Info.i", "ck", "ck.unrestricted", "lk", "uk",
                        "reason.interim", "kMax", "method", "bindingFutility", "cNotBelowFixedc")
        if(any(names(object) %in% valid.args == FALSE)){
            stop("Unknown argument(s) \"",paste(names(object)[names(object) %in% valid.args == FALSE], collapse = "\", \""),"\".\n")
        }
        if(any(valid.args %in% names(object) == FALSE)){
            stop("Missing argument(s) \"",paste(valid.args[valid.args %in% names(object) == FALSE], collapse = "\", \""),"\".\n")
        }

        for(iArg in 1:length(valid.args)){
            assign(valid.args[iArg], value = object[[valid.args[iArg]]], envir = environment())
        }
        alphaSpent <- rep(NA, kMax)
    }else{
        stop("Argument \'object\' should either be \n",
             " a delayedGSD object (output of calcBoundaries functions) \n",
             " a list containing the arguments of the function finalPvalue. \n")
    }
    if(kMax==1){
        stop("kMax should be strictly greater than 1. \n")
    }

    if(missing(continuity.correction)){
        continuity.correction <- options$continuity.correction
    }

    ## ** grid of values
    grid <- NULL
    for(k in 1:(kMax-1)){
        grid <- rbind(grid,
                      data.frame(name = paste0("fut.",k,".",1:sum(seq.futility<ck.unrestricted[k])),
                                 z = seq.futility[seq.futility<ck.unrestricted[k]],
                                 k = k,
                                 alphaSpent = NA)
                      )
        if(cNotBelowFixedc){
            grid <- rbind(grid,
                          data.frame(name = paste0("fut.",k,".Bl"), z = ck.unrestricted[k], k = k, alphaSpent = NA),
                          data.frame(name = paste0("fut.",k,".Bu"), z = 0.1*ck.unrestricted[k]+0.9*ck[k], k = k, alphaSpent = NA))
        }else{
            grid <- rbind(grid,
                          data.frame(name = paste0("fut.",k,".B"), z = ck[k]-1e-5, k = k, alphaSpent = NA))
        }
    }

    grid <- rbind(grid,
                  data.frame(name = paste0("fut.",kMax,".",1:sum(seq.kMax<ck[kMax])),
                             z = seq.kMax[seq.kMax<ck[kMax]],
                             k = kMax,
                             alphaSpent = NA),
                  data.frame(name = paste0("eff.",kMax,".B"),
                             z = ck[kMax],
                             k = kMax,
                             alphaSpent = alphaSpent[kMax]),
                  data.frame(name = paste0("eff.",kMax,".",1:sum(seq.kMax>ck[kMax])),
                             z = seq.kMax[seq.kMax>ck[kMax]],
                             k = kMax,
                             alphaSpent = NA)
                  )
    for(k in (kMax-1):1){
        grid <- rbind(grid,
                      data.frame(name = paste0("eff.",k,".B"),
                                 z = ck[k]+1e-5,
                                 k = k,
                                 alphaSpent = alphaSpent[k]),
                      data.frame(name = paste0("eff.",k,".",1:sum(seq.efficacy>ck[k])),
                                 z = seq.efficacy[seq.efficacy>ck[k]],
                                 k = k,
                                 alphaSpent = NA)
                      )
        
    }
    n.grid <- NROW(grid)

    ## ** evaluate p-value over the domain
    calcP <- function(z, k){
        out <- do.call(FCT.p_value, list(Info.d = Info.d[1:k],
                                         Info.i = Info.i[1:min(k,kMax-1)],
                                        #ck = ck[1:min(k,kMax-1)], ck.unrestricted = ck.unrestricted[1:min(k,kMax-1)],
                                         ck = ck[1:min(k,kMax)],
                                         ck.unrestricted = ck.unrestricted[1:min(k,kMax)],
                                         lk = lk[1:min(k,kMax-1)],
                                         uk = uk[1:min(k,kMax-1)],
                                         kMax = kMax,
                                         delta = 0, 
                                         estimate = z / sqrt(Info.d[k]),
                                         reason.interim = reason.interim[1:k],
                                         method = method,
                                         bindingFutility = bindingFutility, 
                                         cNotBelowFixedc = cNotBelowFixedc,
                                         continuity.correction = continuity.correction)
                       )
        return(data.frame(z = z, k = k, p.value = out, continuity.correction = continuity.correction))
    }
    out <- do.call(rbind,lapply(1:n.grid, function(iG){ ## iG <- 1
        cbind(calcP(z = grid[iG,"z"], k = grid[iG,"k"]), alphaSpent = grid[iG,"alphaSpent"])
    }))

    ## ** graphical display
    if(plot){
        .plotFinalPvalue(out, cNotBelowFixedc = cNotBelowFixedc,
                         col = col, title = title, xlim = xlim, pch = pch, cex = cex, lwd = lwd, digits = digits)
    }

    ## ** export
    if(plot){
        return(invisible(out))
    }else{
        return(out)
    }
}

## * plotFinalPvalue
.plotFinalPvalue <- function(object,
                             cNotBelowFixedc,
                             col,
                             title,
                             xlim,
                             pch,
                             cex,
                             lwd,
                             digits){

    ## ** graphical display
    kMax <- max(object$k)
    if(is.null(xlim)){
        xlim <- c(1-0.1*digits,kMax+0.1*digits)
    }
    n.col <- 1 + 4*(kMax-1)
    if(is.null(col)){
        col <- hcl.colors("Blues", n = n.col+3)
    }else if(length(col)!=n.col){
        stop("Argument color should have length ",n.col,". \n")
    }
    index.change0 <- which(diff(object$k)!=0)
    index.boundary <- c(index.change0[1:(kMax-1)],which.min(abs(object$p.value-0.025)))
    object$p.display <- paste0(round(100*object$p.value, digits = digits),"%")


    plot(object$k,object$z, type = "l", xlab = "stage", ylab = "test statistic", axes = FALSE, xlim = xlim, main = title)
    for(iStage in 1:(kMax-1)){
        if(cNotBelowFixedc){
            points(object$k[index.boundary[iStage]-1],object$z[index.boundary[iStage]-1], pch = pch, cex = cex)
            points(object$k[index.boundary[iStage]],object$z[index.boundary[iStage]], pch = pch, cex = cex)
        }else{
            points(object$k[index.boundary[iStage]],object$z[index.boundary[iStage]], pch = pch, cex = cex)
        }
    }
    points(object$k[index.boundary[kMax]],object$z[index.boundary[kMax]], pch = pch, cex = cex)
    axis(1, at = 0:kMax)
    axis(2)
    for(iStage in 1:(kMax-1)){
        if(cNotBelowFixedc){
            arrows(x0 = object[c(1,index.change0+1)[iStage],"k"], y0 = object[c(1,index.change0+1)[iStage],"z"], 
                   x1 = object[index.change0[iStage]-1,"k"], y1 = object[index.change0[iStage]-1,"z"],
                   col = col[1+2*(iStage-1)], lwd = lwd)
            points(x = object[c(index.change0[iStage]-1, index.change0[iStage]),"k"], y = object[c(index.change0[iStage]-1, index.change0[iStage]),"z"],
                   col = "red", lwd = lwd, type = "l")
        }else{
            arrows(x0 = object[c(1,index.change0+1)[iStage],"k"], y0 = object[c(1,index.change0+1)[iStage],"z"], 
                   x1 = object[index.change0[iStage],"k"], y1 = object[index.change0[iStage],"z"],
                   col = col[1+2*(iStage-1)], lwd = lwd)
        }
        arrows(x0 = object[index.change0[iStage],"k"], y0 = object[index.change0[iStage],"z"],
               x1 = object[index.change0[iStage]+1,"k"], y1 = object[index.change0[iStage]+1,"z"],
               col = col[2+2*(iStage-1)], lwd = lwd)
        iSeq <- seq(c(1,index.change0+1)[iStage],index.change0[iStage]-1,by=1)
        text(x = object[iSeq,"k"], y = object[iSeq,"z"],label = object[iSeq,"p.display"], pos = 2)
        text(x = object[index.change0[iStage],"k"], y = 0.35*object[max(iSeq),"z"]+0.65*object[index.change0[iStage],"z"],
             label = object[index.change0[iStage],"p.display"], pos = 2, col = ifelse(cNotBelowFixedc,"red","black"))
    }
    arrows(x0 = object[index.change0[kMax-1]+1,"k"], y0 = object[index.change0[kMax-1]+1,"z"],
           x1 = object[index.change0[kMax],"k"], y1 = object[index.change0[kMax],"z"],
           col = col[1+2*(kMax-1)], lwd = lwd)
    text(x = object[(index.boundary[kMax-1]+1):index.change0[kMax],"k"], y = object[(index.boundary[kMax-1]+1):index.change0[kMax],"z"],
         label = object[(index.boundary[kMax-1]+1):index.change0[kMax],"p.display"], pos = 4)

    for(iStage in 1:(kMax-1)){ ## iStage <- 1
        arrows(x0 = object[index.change0[kMax+iStage-1],"k"], y0 = object[index.change0[kMax+iStage-1],"z"],
               x1 = object[index.change0[kMax+iStage-1]+1,"k"], y1 = object[index.change0[kMax+iStage-1]+1,"z"],
               col = col[1+2*(kMax+iStage-1)], lwd = lwd)
        arrows(x0 = object[index.change0[kMax+iStage-1]+1,"k"], y0 = object[index.change0[kMax+iStage-1]+1,"z"],
               x1 = object[c(index.change0,NROW(object))[kMax+iStage],"k"], y1 = object[c(index.change0,NROW(object))[kMax+iStage],"z"],
               col = col[2+2*(kMax+iStage-1)], lwd = lwd)
        text(x = object[index.change0[kMax+iStage-1]+1,"k"], y = object[index.change0[kMax+iStage-1]+1,"z"]*1.05,
             label = object[index.change0[kMax+iStage-1]+1,"p.display"], pos = 2)
        iSeq <- seq(index.change0[kMax+iStage-1]+2,c(index.change0,NROW(object))[kMax+iStage],by=1)
        text(x = object[iSeq,"k"], y = object[iSeq,"z"], label = object[iSeq,"p.display"], pos = 2)
    }

    return(invisible(NULL))
}


##----------------------------------------------------------------------
### plotFinalPvalue.R ends here
