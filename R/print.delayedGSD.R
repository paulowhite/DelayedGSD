## * print.delayedGSD
#' @export
print.delayedGSD <- function(x, planned = TRUE, digits = 3, space = " ", abreviated = TRUE, ...){

    ## ** extract information from object
    call <- x$call
    kMax <- x$kMax
    k <- x$stage$k
    type.k <- x$stage$type
    InflationFactor <- x$InflationFactor

    test.planning <- (type.k=="planning") || identical(planned,"only")
    ls.info <- getInformation(x, planned = planned)
    Info.max <- ls.info$Info.max
    uk  <- ls.info$uk
    lk  <- ls.info$lk
    ck  <- ls.info$ck
    
    if(!test.planning){
        iLMM <- x$lmm[[utils::tail(ls.info$index.lmm,1)]]

        statistic <- unname(ls.info$delta[,"statistic"])
        statistic.interim <- c(statistic[1:k], rep(NA,kMax-k))
        statistic.interim[kMax] <- NA ## at final analysis no interim only decision
        statistic.decision <- rep(NA,kMax)
        if(type.k=="decision"){
            statistic.decision[k] <- statistic[k+1]
        }else if(type.k=="final"){
            statistic.decision[k] <- statistic[k]
        }

        ## remove critical boundaries when continuing at interim
        ck[which(x$conclusion["interim",]=="continue")] <- NA
    }

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
    if(test.planning){

        df.print <- data.frame(1:kMax, lk, uk, c(ck,uk[kMax]))
        if(abreviated){
            colnames(df.print) <- c("Stage","F-bound","E-bound","C-bound")
        }else{
            colnames(df.print) <- c("Stage","Futility boundary","Efficacy boundary","Critical boundary")
        }

        cat("\n * Planned boundaries: \n")
        print(df.print, row.names = FALSE, digits = digits)

    }else{    

        ## Round values
        if(is.null(digits)){
            statistic.interim.print <- as.character(statistic.interim)
            statistic.decision.print <- as.character(statistic.decision)
            lk.print <- c(as.character(lk[1:(kMax-1)]),as.character(NA)) ## no interim at the last stage, only decision
            uk.print <- c(as.character(uk[1:(kMax-1)]),as.character(NA))
            ck.print <- as.character(c(ck,uk[kMax]))
        }else{
            statistic.interim.print <- as.character(round(statistic.interim,digits))
            statistic.decision.print <- as.character(round(statistic.decision,digits))
            lk.print <- c(as.character(round(lk[1:(kMax-1)], digits)),as.character(NA)) ## no interim at the last stage, only decision
            uk.print <- c(as.character(round(uk[1:(kMax-1)], digits)),as.character(NA))
            ck.print <- c(as.character(round(c(ck,uk[kMax]),digits)))
        }
        statistic.interim.print[is.na(statistic.interim)] <- ""
        statistic.decision.print[is.na(statistic.decision)] <- ""
        lk.print[is.na(lk.print)] <- ""
        uk.print[is.na(uk.print)] <- ""
        ck.print[is.na(ck.print)] <- ""

        ## Abreviate actions
        colname.interim <- c("Futility boundary","Efficacy boundary","Statistic")
        colname.decision <- c("Critical boundary","Statistic")
        
        conclusion.interim.print <- .reformatInterimTxt(x$conclusion, kMax = kMax, abreviated = abreviated)
        if(!is.null(conclusion.interim.print)){
            if(abreviated){
                colname.interim <- c(colname.interim,"")
            }else{
                colname.interim <- c(colname.interim,"Action")
            }
        }
        conclusion.decision.print <- .reformatDecisionTxt(x$conclusion, kMax = kMax, abreviated = abreviated)
        if(!is.null(conclusion.decision.print)){
            if(abreviated){
                colname.decision <- c(colname.decision,"")
            }else{
                colname.decision <- c(colname.decision,"Action")
            }
        }

        ## assemble
        if(is.null(space)){
            ls.df <- list(1:kMax,
                          lk.print, uk.print, statistic.interim.print, conclusion.interim.print,
                          ck.print, statistic.decision.print, conclusion.decision.print)
            df.print <- data.frame(ls.df[sapply(ls.df,length)>0])
            colnames(df.print) <- c("Stage",colname.interim,colname.decision)
        }else{
            ls.df <- list(1:kMax, space,
                          lk.print, uk.print, statistic.interim.print, conclusion.interim.print, space,
                          ck.print, statistic.decision.print, conclusion.decision.print)
            df.print <- data.frame(ls.df[sapply(ls.df,length)>0])
            colnames(df.print) <- c("Stage","",colname.interim,"",colname.decision)
        }

        ## display
        cat("\n * Boundaries and observed statistics \n")
        df.print2 <- rbind(names(df.print),df.print)
        colnames(df.print2) <- rep("", NCOL(df.print))
        colnames(df.print2)[names(df.print)=="Efficacy boundary"] <- "Interim"
        colnames(df.print2)[names(df.print)=="Critical boundary"] <- "Decision"

        if(abreviated){
            df.print2[1,names(df.print)=="Statistic"] <- "Stat"
            df.print2[1,names(df.print)=="Futility boundary"] <- "F-bound"
            df.print2[1,names(df.print)=="Efficacy boundary"] <- "E-bound"
            df.print2[1,names(df.print)=="Critical boundary"] <- "C-bound"
        }
        print(df.print2, row.names = FALSE)
    }
  
    cat("\n")

    ## ** Information
    if(test.planning){
        Info.i <- c(x$planned$Info.i[1:(kMax-1)],NA)
        Info.d <- c(x$planned$Info.d,x$planned$Info.i[kMax])
        InfoR.i <- Info.i/x$Info.max
        InfoR.d <- Info.d/x$Info.max
        if(!is.null(x$n.obs)){
            n.obs <- x$n.obs*InfoR.d
        }else{
            n.obs <- NULL
        }
        cat(" * Planned information: \n",sep="")
    }else{
        Info.i <- c(ls.info$Info.i[1:(kMax-1)],NA)
        Info.d <- c(ls.info$Info.d,ls.info$Info.i[kMax])
        InfoR.i <- Info.i/x$Info.max
        InfoR.d <- Info.d/x$Info.max
        n.obs <- unlist(lapply(x$lmm, function(iLMM){
            if(!is.null(iLMM)){iLMM$n["decision"]}else{NA}
        }))
        index.lastNNA <- utils::tail(which(!is.na(n.obs)),1)
        index.NA <- which(is.na(n.obs))
        n.obs[index.NA] <- round(n.obs[index.lastNNA] * InfoR.d[index.NA]/InfoR.d[index.lastNNA],digits)
        n.obs <- as.character(n.obs)
        n.obs[intersect(index.NA,which(is.na(InfoR.d)))] <- ""
        cat(" * Observed and planned information: \n",sep="")
        ## cat(" * Observed and predicted information: \n",sep="")
    }
    if(!is.null(digits)){
        Info.i <- round(Info.i, digits)
        Info.d <- round(Info.d, digits)
        InfoR.i <- round(InfoR.i, digits)
        InfoR.d <- round(InfoR.d, digits)
    }
    df.printInfo <- data.frame(Stage = 1:kMax,
                               Interim =  as.character(Info.i),
                               "(%)" = as.character(InfoR.i),
                               "Decision" = as.character(Info.d),
                               "(%)" = as.character(InfoR.d),
                               check.names = FALSE)
    df.printInfo[is.na(df.printInfo)] <- ""
    if(!is.null(n.obs)){
        df.printInfo  <- cbind(df.printInfo,"n"=n.obs)
    }
    if(!is.null(space)){
        df.printInfo <- cbind(df.printInfo[,1,drop=FALSE]," "=space,df.printInfo[,2:3]," "=space,df.printInfo[,-(1:3)])
        names(df.printInfo)[c(2,5)] <- ""
    }
        

    print(df.printInfo, row.names = FALSE)
    cat("\n")

    ## ** Estimates
    if(!test.planning){
        cat(" * Current ML-estimate of the treatment effect \n")
        print(c("estimate" = iLMM$estimate, "se" = iLMM$se, "statistic" = iLMM$statistic, "df" = iLMM$df, "p.value" = iLMM$p.value), digits = digits)
        cat("\n")
    }  
  
    ## ** Sample size
    if(test.planning){
        cat(" * Inflation factor: ",InflationFactor," \n",sep="")
    }else{
        cat(" * Number of clusters in the study: ",iLMM$n["interim"]," with at least one outcome value \n ",
            "                                   ",iLMM$n["interim.cc"]," with complete data \n", sep = "")
    }    
    return(invisible(NULL))
}

## * reformat text
## ** interim
.reformatInterimTxt <- function(conclusion, kMax, abreviated){

    if(all(is.na(conclusion["interim",]))){return(NULL)}
    
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
            }else if(conclusion["reason.interim",iK]=="Imax-reached"){
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

    if(all(is.na(conclusion["decision",]))){return(NULL)}

    out <- sapply(1:kMax,function(iK){
        if(is.na(conclusion["decision",iK])){
            return("")
        }else if(conclusion["decision",iK]=="Efficacy"){
            if(abreviated){
                return("E")
            }else{
                return("conclude efficacy")
            }
        }else if(conclusion["decision",iK]=="Futility"){
            if(abreviated){
                return("F")
            }else{
                return("conclude futility")
            }
        }
    })

    return(out)
}
