### AnalyzeData.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 25 2020 (09:23) 
## Version: 
## Last-Updated: mar 11 2021 (12:06) 
##           By: Brice Ozenne
##     Update #: 78
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

#' @param d dataset
#' @param ddf [character] method used to compute the degrees of freedom of the Wald statistic.
#' Can be \code{satterthwaite} or \code{nlme}.
AnalyzeData <- function(d, ddf = "nlme"){

    require(nlme)

    ## ** normalize arguments
    ddf <- match.arg(ddf, choices = c("nlme","satterthwaite"))
    long <- wide2long(d)
    ## head(long[order(long$id),],n=10)
    ## summary(long)

    ## ** fit gls model
    ctrl <- glsControl(opt='optim')
    
    m <- gls(X~baseline*time.factor + Z*time.factor,
             data=long,
             correlation=corSymm(form=~1|id), 
             varIdent(form=~1|time),method="REML",
             na.action = na.exclude)

    ## ** store estimate from the gls model    
    out <- list(d.long=long,
                fit=m,
                pvalue=NA,
                estimate=coef(m)["Z1"],
                se=sqrt(m$varBeta["Z1","Z1"])
                )

    ## computation of the p-value  
    if(ddf=="nlme"){ ## with weird estimator of the degree of freedom from nlme::gls 
        iStat <- abs(coef(m)["Z1"]/sqrt(m$varBeta["Z1","Z1"]))
        iDF <- m$dims$N - m$dims$p
        out$pvalue <- 2*(1-pt(iStat, df = iDF))
        ## summary(m)$tTable["Z1","p-value"]
    }else{ ## or using satterthwaite approximation
        require(emmeans)
        groupTest <- emmeans::emmeans(m, specs = ~Z|time.factor, data = long[!is.na(long$X),,drop=FALSE])
        e.satterthwaite <- summary(pairs(groupTest, reverse = TRUE), by = NULL, infer = TRUE, adjust = "none")
        index <- which(e.satterthwaite$time.factor==levels(long$time.factor)[1])

        if(abs(out$estimate-coef(m)["Z1"])>1e-10){
            warning("Something went wrong when extracting the coefficient from emmeans \n",
                    "discrepancy estimate emmeans vs. gls: ",out$estimate-coef(m)["Z1"])
        }
        if(abs(1/e.satterthwaite$SE[index]^2-out$Info)>1e-10){
            warning("Something went wrong when estimating the information in emeans \n",
                    "discrepancy se emmeans vs. gls: ",1/e.satterthwaite$SE[index]^2-out$Info)
        }
        out$p.value <- e.satterthwaite$p.value[index]

    }


    ## ** Estimate the information
    ## current information and information had we had no missing values
    res.info <- attr(getInformation(m, name.coef = "Z1", type = "prediction", method.prediction = "inflation",
                               data = long, details = TRUE), "details")


    if(any(abs(res.info$interim$vcov - vcov(m))>1e-10)){
        warning("Something went wrong when extracting the variance covariance matrix from gls \n",
                "largest discrepancy between getInformation and gls: ",max(abs(res.info$interim$vcov - vcov(m))))
     
    }
    out$Info <- res.info$interim$information
    out$Info.fullData <- res.info$decision$information
    out$n.cc <- res.info$n["interim.cc"]
    out$n <- res.info$n["interim"]
    out$n.fullData <- res.info$n["decision"]
    
    ## ** Export    
    return(out)
}

wide2long <- function(d, rm.na = FALSE, id.na = NULL){

    ## get the number of timepoints
    N.fw <- length(grep("t",names(d)))
    
    ## rename column with baseline value
    names(d)[names(d)=="X1"] <- "baseline"

    ## convert to long format
    long <- reshape(d, direction='long', varying=paste0("X",2:N.fw), idvar='id', v.names='X')
    long$id <- factor(long$id)
    long$Z <- factor(long$Z)
    long$time.factor <- factor(long$time)
    long$time.factor <- relevel(long$time.factor,ref=N.fw-1)## relvel to make the treatment coefficient have the interpretation of threatment effect at end of follow-up

    ## remove lines corresponding to dropout
    if(rm.na){
        long <- long[!is.na(long$X),]
    }

    ## add missing values for not yet observed ids
    if(!is.null(id.na)){
        long[long$id %in% id.na, c("baseline","t1","t2","t3","X")] <- NA
    }
    
    ## export
    return(long)
}
##----------------------------------------------------------------------
### AnalyzeData.R ends here
