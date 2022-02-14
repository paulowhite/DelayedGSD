## * print.delayedGSD (documentation)
#' @title Status of a Group Sequential Design with Delayed Endpoints
#' @description Display boundaries and estimated treatment effect up the current stage.
#' 
#' @param x object of class \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param planned [logical] should the planned boundaries be displayed along with the boundaries that have been computed based on the data.
#' Can also be \code{"only"} to only display planned boundaries.
#' @param digits [interger,>0] integer indicating the number of decimal places to be used.
#' @param space [character] column separator.
#' @param abreviated [logical] abreviate some column names and action that should be taken
#' (e.g. C for continue, C-Info for continue due to decreasing information, S-F for stop for futility, and S-Imax for stop because maximum information has been reached).
#' @param ... not used, for compatibility with the generic method.
#' 

## * print.delayedGSD
#' @export
print.delayedGSD <- function(x, planned = TRUE, digits = 5, space = " ", abreviated = TRUE, ...){

    ## ** extract information from object
    call <- x$call
    kMax <- x$kMax
    k <- x$stage$k
    type.k <- x$stage$type
    
    test.planning <- (type.k=="planning") || identical(planned,"only")

    ## ** Welcome message
    if(test.planning){

        cat("	    Planning for a GSD with repeated measurements \n \n")
        cat(" * Call \n")
        print(call)

    }else{

        if(type.k == "final"){
            cat("	    GSD with repeated measurements at the final analysis (stage ",k,") \n", sep = "")
        }else if(type.k == "decision"){
            cat("	    GSD with repeated measurements at the decision analysis of stage ",k," \n", sep = "")
        }else if(type.k == "interim"){
            cat("	    GSD with repeated measurements at the interim analysis of stage ",k," \n", sep = "")
        }

    }

    ## ** Boundaries
    df.printBound <- stats::coef(x, type = "boundary", planned = planned)
    df.printBound$Fbound <- c(round(df.printBound$Fbound[-NROW(df.printBound)], digits),NA)
    df.printBound$Ebound <- c(round(df.printBound$Ebound[-NROW(df.printBound)], digits),NA)
    df.printBound$Cbound <- round(df.printBound$Cbound, digits)

    if(test.planning){

        if(abreviated){
            colnames(df.printBound) <- c("Stage","F-bound","E-bound","C-bound")
        }else{
            colnames(df.printBound) <- c("Stage","Futility boundary","Efficacy boundary","Critical boundary")
        }
        df.printBound[is.na(df.printBound)] <- ""
        df.printBound[["alpha-spent"]] <- round(x$planned$alphaSpent, digits)
        df.printBound[["beta-spent"]] <- round(x$planned$betaSpent, digits)
        cat("\n * Planned boundaries: \n")
        print(df.printBound, row.names = FALSE, quotes="")

    }else{
        df.printBound$statistic.interim <- round(df.printBound$statistic.interim, digits)
        df.printBound$statistic.decision <- round(df.printBound$statistic.decision, digits)
        if(planned==FALSE){
            df.printBound[["alpha-spent"]] <- NA
            df.printBound[["beta-spent"]] <- NA
            df.printBound[["alpha-spent"]][1:x$stage$k] <- round(x$alphaSpent[1:x$stage$k], digits+2)
            df.printBound[["beta-spent"]][1:x$stage$k] <- round(x$betaSpent[1:x$stage$k], digits+2)
        }else{
            df.printBound[["alpha-spent"]] <- round(x$alphaSpent, digits+2)
            df.printBound[["beta-spent"]] <- round(x$betaSpent, digits+2)
        }
        df.printBound[is.na(df.printBound)] <- ""
        conclusion.interim.print <- .reformatInterimTxt(x$conclusion, kMax = kMax, abreviated = abreviated)
        conclusion.decision.print <- .reformatDecisionTxt(x$conclusion, kMax = kMax, abreviated = abreviated)

        df.printBound2 <- data.frame("Stage" = df.printBound$Stage,
                                     "Futility boundary" = df.printBound$Fbound,
                                     "Efficacy boundary" = df.printBound$Ebound,
                                     "Statistic" = df.printBound$statistic.interim,
                                     "Action" = conclusion.interim.print,
                                     "Critical boundary" = df.printBound$Cbound,
                                     "Statistic" = df.printBound$statistic.decision,
                                     "Action" = conclusion.decision.print,
                                     "alpha" = df.printBound[["alpha-spent"]],
                                     "beta" = df.printBound[["beta-spent"]],
                                     check.names = FALSE)

        if(abreviated){
            names(df.printBound2) <- c("Stage", "F-bound", "E-bound", "Stat", "", "C-bound", "Stat", "","alpha","beta")
        }

        ## display
        cat("\n * Boundaries and observed statistics \n")
        df.printBound3 <- rbind(names(df.printBound2),df.printBound2)
        df.printBound3[1,1] <- ""
        names(df.printBound3) <- rep("", NCOL(df.printBound3))
        if(abreviated){
            names(df.printBound3)[3] <- "Interim"
            names(df.printBound3)[6] <- "Decision"
            names(df.printBound3)[9] <- "Spent"
        }else{
            names(df.printBound3)[3] <- "Interim analysis"
            names(df.printBound3)[6] <- "Decision analysis"
            names(df.printBound3)[9] <- "Error spent"
        }
        print(df.printBound3, row.names = FALSE, quote = FALSE)
    }
    cat("\n")

    ## ** Information
    df.printInfo <- stats::coef(x, type = "information", planned = planned)
    df.printInfo$Interim[length(df.printInfo$Interim)] <- NA
    df.printInfo$Interim.pc[length(df.printInfo$Interim.pc)] <- NA
    if(!is.null(digits)){        
        df.printInfo$Interim <- round(df.printInfo$Interim, digits)
        df.printInfo$Decision <- round(df.printInfo$Decision, digits)
        df.printInfo$Interim.pc <- round(df.printInfo$Interim.pc, digits)
        df.printInfo$Decision.pc <- round(df.printInfo$Decision.pc, digits)
    }

    if(test.planning){
        names(df.printInfo) <- c("Stage", "Interim", "(%)", "Decision", "(%)")
        cat(" * Planned information: \n",sep="")
    }else{
        names(df.printInfo) <- c("Stage", "Interim", "(%)", "Decision", "(%)", "n")
        if(planned){
            cat(" * Observed and planned information: \n",sep="")
        }else{
            cat(" * Observed and predicted information: \n",sep="")
        }
    }
    df.printInfo[is.na(df.printInfo)] <- ""
    if(!is.null(space)){
        df.printInfo <- cbind(df.printInfo[,1,drop=FALSE]," "=space,df.printInfo[,2:3]," "=space,df.printInfo[,-(1:3)])
        names(df.printInfo)[c(2,5)] <- ""
    }
    print(df.printInfo, row.names = FALSE)
    cat("\n")

    ## ** Estimates
    if(!test.planning){
        x.CI <- stats::confint(x, method = c("ML","corrected ML"))
        x.CI$estimate <- round(x.CI$estimate, digits)
        x.CI$se <- round(x.CI$se, digits)
        x.CI$lower <- round(x.CI$lower, digits)
        x.CI$upper <- round(x.CI$upper, digits)
        x.CI$statistic <- round(x.CI$statistic, digits)
        x.CI$df <- round(x.CI$df, digits)
        x.CI$p.value <- format.pval(x.CI$p.value, digits)
        x.CI$stage <- paste0(x.CI$stage," (",x.CI$type,")")
        x.CI$type <- NULL
        if(is.null(x$correction)){
            cat(" * Current ML-estimate of the treatment effect (",x.CI$coef[1],") \n",sep="")
            print(x.CI[,c("estimate","se","lower","upper","statistic","df","p.value")], row.names = FALSE)
            cat("\n")
        }else{
            cat(" * Estimate of the treatment effect (",x.CI$coef[1],") \n",sep="")
            x.CI[is.na(x.CI)] <- ""
            
            print(x.CI[,c("method","estimate","lower","upper","p.value")], quote = FALSE, row.names = FALSE)
            cat("\n")
        }
    }  
  
    ## ** Sample size
    if(test.planning){
        InflationFactor <- x$planned$InflationFactor
        cat(" * Inflation factor: ",InflationFactor," \n",sep="")
    }else{
        iLMM <- x$lmm[[k+(type.k=="decision")]]
        cat(" * Number of clusters in the study: ",iLMM$n["interim"]," with at least one outcome value \n ",
            "                                   ",iLMM$n["interim.cc"]," with complete data \n", sep = "")
    }    
    return(invisible(NULL))
}

## * reformat text
## ** interim
.reformatInterimTxt <- function(conclusion, kMax, abreviated){

    if(all(is.na(conclusion["interim",]))){return(rep("",kMax))}
    
    out <- sapply(1:kMax,function(iK){
        if(is.na(conclusion["interim",iK])){
            return("")
        }else if(conclusion["interim",iK]=="continue"){
            if(conclusion["reason.interim",iK]=="no boundary crossed"){
                if(abreviated){
                    return("C")
                }else{
                    return("continue")
                }
            }else if(conclusion["reason.interim",iK]=="decreasing information"){
                if(abreviated){
                    return("C-Info")
                }else{
                    return("continue (decreasing information)")
                }
            }
        }else if(conclusion["interim",iK]=="stop"){
            if(conclusion["reason.interim",iK]=="efficacy"){
                if(abreviated){
                    return("S-E")
                }else{
                    return("stop for efficacy")
                }
            }else if(conclusion["reason.interim",iK]=="futility"){
                if(abreviated){
                    return("S-F")
                }else{
                    return("stop for futility")
                }
            }else if(conclusion["reason.interim",iK]=="Imax reached"){
                if(abreviated){
                    return("S-Imax")
                }else{
                    return("stop (maximum information reached)")
                }
            }
        }
    })
    return(out)
}

## ** decision
.reformatDecisionTxt <- function(conclusion, kMax, abreviated){

    out <- rep("",kMax)

    if(all(is.na(conclusion["decision",]))){return(out)}

    out <- sapply(1:kMax,function(iK){
        if(is.na(conclusion["decision",iK])){
            return("")
        }else if(conclusion["decision",iK]=="efficacy"){
            if(abreviated){
                return("E")
            }else{
                return("conclude efficacy")
            }
            return(out)
        }else if(conclusion["decision",iK]=="futility"){
            if(abreviated){
                return("F")
            }else{
                return("conclude futility")
            }
        }
    })

    return(out)
}
