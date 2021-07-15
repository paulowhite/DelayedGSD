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
        statistic.decision <- rep(NA,kMax)
        if(type.k=="decision"){
            statistic.decision[k] <- statistic[k+1]
        }else if(type.k=="final"){
            statistic.decision[k] <- statistic[k]
        }

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

    ## ** Boundaries and progress of the trial
    if(test.planning){

        df.print <- data.frame(1:kMax, lk, uk, c(ck,uk[kMax]))
        if(abreviated){
            colnames(df.print) <- c("Stage","F-bound","E-bound","C-bound")
        }else{
            colnames(df.print) <- c("Stage","Futility boundary","Efficacy boundary","Critical boundary")
        }

        cat("\n * Boundaries \n")
        print(df.print, row.names = FALSE, digits = digits)

    }else{    

        ## convert to character interim
        if(is.null(digits)){
            statistic.interim.print <- as.character(statistic.interim)
            lk.print <- lk
            uk.print <- uk
        }else{
            statistic.interim.print <- as.character(round(statistic.interim,digits))
            lk.print <- round(lk, digits)
            uk.print <- round(uk, digits)
        }
        statistic.interim.print[is.na(statistic.interim)] <- ""
        colname.interim <- c("Futility boundary","Efficacy boundary","Statistic")
        
        if(any(!is.na(x$conclusion["interim",]))){
            conclusion.interim.print <- sapply(1:kMax,function(iK){
                if(is.na(x$conclusion["interim",iK])){
                    return("")
                }else if(x$conclusion["interim",iK]=="Efficacy"){
                    if(abreviated){
                        return("E")
                    }else{
                        return("conclude efficacy")
                    }
                }else if(x$conclusion["interim",iK]=="Futility"){
                    if(abreviated){
                        return("F")
                    }else{
                        return("conclude futility")
                    }
                }else if(x$conclusion["interim",iK]=="continue"){
                    if(x$conclusion["reason.interim",iK]=="no boundary crossed"){
                        if(abreviated){
                            return("C")
                        }else{
                            return("continue")
                        }
                    }else if(x$conclusion["reason.interim",iK]=="decreasing information"){
                        if(abreviated){
                            return("C-Info")
                        }else{
                            return("continue (decreasing information)")
                        }
                    }
                }else if(x$conclusion["interim",iK]=="stop"){
                    if(x$conclusion["reason.interim",iK]=="efficacy"){
                        if(abreviated){
                            return("S-E")
                        }else{
                            return("stop for efficacy")
                        }
                    }else if(x$conclusion["reason.interim",iK]=="futility"){
                        if(abreviated){
                            return("S-F")
                        }else{
                            return("stop for futility")
                        }
                    }else if(x$conclusion["reason.interim",iK]=="Imax-reached"){
                        if(abreviated){
                            return("S-Imax")
                        }else{
                            return("stop (maximum information reached)")
                        }
                    }
                }
            })
            colname.interim <- c(colname.interim,"")
        }else{
            conclusion.interim.print <- NULL
        }

        ## convert to character decision
        if(is.null(digits)){
            statistic.decision.print <- as.character(statistic.decision)
            ck.print <- ck
        }else{
            statistic.decision.print <- as.character(round(statistic.decision,digits))
            ck.print <- round(ck,digits)
        }
        statistic.decision.print[is.na(statistic.decision)] <- ""
        colname.decision <- c("Critical boundary","Statistic")
        if(any(!is.na(x$conclusion["decision",]))){
            conclusion.decision.print <- x$conclusion["decision",]
            if(abreviated){
                conclusion.decision.print[conclusion.decision.print=="Efficacy"] <- "E"
                conclusion.decision.print[conclusion.decision.print=="Futility"] <- "F"
            }else{
                conclusion.decision.print[conclusion.decision.print=="Efficacy"] <- "conclude efficacy"
                conclusion.decision.print[conclusion.decision.print=="Futility"] <- "conclude futility"
            }
            conclusion.decision.print[is.na(x$conclusion["decision",])] <- ""
            colname.decision <- c(colname.decision,"")
        }else{
            conclusion.decision.print <- NULL
        }

        ## assemble
        if(is.null(space)){
            ls.df <- list(1:kMax,
                          lk.print, uk.print, statistic.interim.print, conclusion.interim.print,
                          c(ck.print,uk.print[kMax]), statistic.decision.print, conclusion.decision.print)
            df.print <- data.frame(ls.df[sapply(ls.df,length)>0])
            colnames(df.print) <- c("Stage",colname.interim,colname.decision)
        }else{
            ls.df <- list(1:kMax, space,
                          lk.print, uk.print, statistic.interim.print, conclusion.interim.print, space,
                          c(ck.print,uk.print[kMax]), statistic.decision.print, conclusion.decision.print)
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

    if(test.planning){
        cat(" * Planned maximum information: ",Info.max," \n \n",sep="")
        cat(" * Inflation factor: ",InflationFactor," \n",sep="")
    }else{
        cat(" * Number of clusters: ",iLMM$n["interim"]," with at least one outcome value \n ",
            "                      ",iLMM$n["interim.cc"]," with complete data \n", sep = "")
        cat("\n")
        cat(" * Current estimated treatment effect \n")
        print(c("estimate" = iLMM$estimate, "se" = iLMM$se, "statistic" = iLMM$statistic, "df" = iLMM$df, "p.value" = iLMM$p.value))
        cat("\n")
        cat(" * Current estimated information: ",iLMM$getInformation["interim"]," \n",sep="")
    }  
    cat("\n")
  
    return(invisible(NULL))
}
