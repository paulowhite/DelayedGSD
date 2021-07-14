## * print.delayedGSD
#' @export
print.delayedGSD <- function(x, planned = TRUE, ...){

    ## ** extract information from object
    call <- x$call
    kMax <- x$kMax
    k <- x$stage$k
    decision <- x$stage$decision
    InflationFactor <- x$InflationFactor

    test.planning <- (x$stage$k==0) || identical(planned,"only")
    ls.info <- getInformation(x, planned = planned)
    Info.max <- ls.info$Info.max
    uk  <- ls.info$uk
    lk  <- ls.info$lk
    ck  <- ls.info$ck

    if(!test.planning){
        delta  <- ls.info$delta
        iLMM <- x$lmm[[utils::tail(ls.info$index.lmm,1)]]
    }

    ## ** display
    if(test.planning){
        cat("	    Planning for a GSD with repeated measurements \n \n")
        cat(" * Call \n")
        print(call)

    }else{

        if(k == kMax){
            cat("	    GSD with repeated measurements at the final analysis (stage ",k,") \n", sep = "")
        }else if(decision>0){
            cat("	    GSD with repeated measurements at the decision analysis of stage ",k," \n", sep = "")
        }else{
            cat("	    GSD with repeated measurements at the interim analysis of stage ",k," \n", sep = "")
        }
    }
    cat("\n")

    bnds <- cbind(Stage = c(1:kMax), Lower = lk[1:kMax], Upper = uk[1:kMax], Decision = c(ck,uk[kMax]))

    cat(" * Boundaries \n")
    print(as.data.frame(bnds), row.names = FALSE)
  
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
