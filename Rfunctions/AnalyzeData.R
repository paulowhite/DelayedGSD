### AnalyzeData.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 25 2020 (09:23) 
## Version: 
## Last-Updated: mar 26 2021 (23:48) 
##           By: Brice Ozenne
##     Update #: 195
#----------------------------------------------------------------------
## 
### Commentary: 
##
## Analyze available data: provide estimate, s.e., Information nd p-value.
##
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * AnalyzeData
#' @param d dataset
#' @param ddf [character] method used to compute the degrees of freedom of the Wald statistic.
#' Can be \code{satterthwaite} or \code{nlme}.
#' @param getinfo [logical] should the information be computed at interim and decision? Otherwise no information is computed.
AnalyzeData <- function(d, ddf = "nlme", getinfo = TRUE){

    require(nlme)

    ## ** normalize arguments
    ddf <- match.arg(ddf, choices = c("nlme","satterthwaite"))
    long <- wide2long(d, rm.na = TRUE) ## remove missing value but not pipeline value from the data.set
    ## head(long[order(long$id),],n=10)
    ## summary(long)

    ## ** fit gls model
    ctrl <- glsControl(opt='optim')
    
    m <- gls(X~baseline*visit + Z*visit,
             data=long,
             correlation=corSymm(form=~visit.num|id), 
             varIdent(form=~1|visit),method="REML",
             na.action = na.exclude)

    ## ** store estimate from the gls model    
    out <- list(d.long=long,
                fit=m,
                estimate = as.double(coef(m)["Z1"]),
                se = as.double(sqrt(m$varBeta["Z1","Z1"])),
                Info = 1/m$varBeta["Z1","Z1"], 
                df = NA,
                p.value = NA
                )

    ## computation of the p-value  
    if(ddf=="nlme"){ ## with weird estimator of the degree of freedom from nlme::gls 
        iStat <- abs(out$estimate)/out$se
        out$df <- as.double(m$dims$N - m$dims$p)
        out$p.value <- as.double(2*(1-pt(iStat, df = out$df)))
        ## summary(m)$tTable["Z1","p-value"]
    }else{ ## or using satterthwaite approximation
        require(emmeans)
        groupTest <- emmeans::emmeans(m, specs = ~Z|visit, data = long[!is.na(long$X),,drop=FALSE])
        e.satterthwaite <- summary(pairs(groupTest, reverse = TRUE), by = NULL, infer = TRUE, adjust = "none")
        index <- which(e.satterthwaite$visit==levels(long$visit)[1])

        if(abs(out$estimate - e.satterthwaite$estimate[index])>1e-10){
            warning("Something went wrong when extracting the coefficient from emmeans \n",
                    "discrepancy estimate emmeans vs. gls: ",out$estimate-e.satterthwaite$estimate[index])
        }
        if(abs(out$se - e.satterthwaite$SE[index])>1e-10){
            warning("Something went wrong when extracting the standard error in emeans \n",
                    "discrepancy se emmeans vs. gls: ",e.satterthwaite$SE[index]-out$se)
        }
        out$p.value <- e.satterthwaite$p.value[index]
    }

    ## ** Estimate the information
    if(getinfo){

        ## current information and information at decision
        out <- c(out,getInformation(m, name.coef = "Z1", data = long, details = TRUE))
        names(out)[names(out)=="info"] <- "getInformation"
        if(any(abs(out$interim$vcov - vcov(m))>1e-10)){
            warning("Something went wrong when extracting the variance covariance matrix from gls \n",
                    "largest discrepancy between getInformation and gls: ",max(abs(res.info$interim$vcov - vcov(m))))
     
        }
    }

    ## ** Export
    class(out) <- "getInformationGLS"
    return(out)
}

## * print.getInformationGLS
print.getInformationGLS <- function(x, ...){
    cat("       Analysis via the gls function \n \n ")
    cat("  Estimated treatment effect: ",x$name.coef," \n",sep = "")
    print(c("estimate" = x$coef, "se" = x$se, "df" = x$df, "p.value" = x$p.value))
    cat("\n  Number of clusters \n")
    print(x$n)
    cat("\n  Estimated information \n")
    print(x$getInformation)
    cat("\n",
        "total     : all patients including who dropped out early and have no observable outcome \n",
        "decision  : patients with at least one observable outcome, including those with not yet observed values \n",
        "interim   : patients with at least one observed outcome \n",
        "interim.cc: patients with at no missing outcome \n")
}

## * wide2long
wide2long <- function(d, rm.na = FALSE, id.na = NULL, Z.id.na = TRUE){
    require(data.table)
    if(data.table::is.data.table(d)){
        dtW <- data.table::copy(d)
    }else{
        dtW <- data.table::as.data.table(d)
    }
    ## get the number of timepoints
    col.time <- grep("^t([[:digit:]]+)$",names(dtW), value = TRUE)
    col.X <- grep("^X([[:digit:]]+)$*",names(dtW), value = TRUE)
    col.missing <- grep("^missing([[:digit:]]+)$*",names(dtW), value = TRUE)
    if(length(col.missing)==0){
        dtW$missing1 <- as.numeric(is.na(dtW$X1))
        dtW$missing2 <- as.numeric(is.na(dtW$X2))
        dtW$missing3 <- as.numeric(is.na(dtW$X3))
        col.missing <- grep("^missing([[:digit:]]+)$*",names(dtW), value = TRUE)
    }
    n.times <- length(grep("t",names(dtW)))

    ## convert to long format
    dtL <- data.table::melt.data.table(dtW, id.vars = c("id","Z","missing1","t1","X1"),
                                       variable.name = "visit",
                                       measure.vars = list("missing" = c("missing2","missing3"),
                                                           "time" = c("t2","t3"),
                                                           "X" = c("X2","X3")))

    dfL <- as.data.frame(dtL)
    names(dfL)[names(dfL)=="X1"] <- "baseline"
    dfL$id <- factor(dfL$id)
    dfL$Z <- factor(dfL$Z)
    dfL$visit.num <- as.numeric(dfL$visit)
    dfL$visit <- relevel(dfL$visit,ref=n.times-1)## relevel to make the treatment coefficient have the interpretation of threatment effect at end of follow-up
    ## add missing values for not yet observed ids
    if(!is.null(id.na)){
        index.notObs <- which(dfL$id %in% id.na)
        dfL[index.notObs, c("missing1","missing")] <- -1
        dfL[index.notObs, c("t1","baseline","time","X")] <- NA
        if(Z.id.na){
            dfL[index.notObs, c("Z")] <- NA
        }
    }
    dfL.save <- dfL

    table(dfL$missing1)
    ## remove lines corresponding to missing data
    if(rm.na){
        dfL <- dfL[(dfL$missing1<=0),] ## at baseline 
        dfL <- dfL[(dfL$missing<=0),] ## or at follow-up
        dfL$id <- droplevels(dfL$id)
        dfL$Z <- droplevels(dfL$Z)
        dfL$visit <- droplevels(dfL$visit)
    }


    ## export
    attr(dfL,"df.allobs") <- dfL.save
    return(dfL)
}
##----------------------------------------------------------------------
### AnalyzeData.R ends here
