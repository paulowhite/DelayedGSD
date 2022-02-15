## * analyzeData (documentation)
#' @title Fit a Linear Mixed Model
#' @description Fit a linear mixed model and extract the information relative to the parameter of interest.
#' 
#' @param d dataset
#' @param ddf [character] method used to compute the degrees of freedom of the Wald statistic.
#' Can be \code{"satterthwaite"} or \code{"nlme"}.
#' @param getinfo [logical] should the information be computed at interim and decision? Otherwise no information is computed.
#' @param trace [logical] should the execution of the function be traced?
#' 
#' @examples
#' ## simulate data in the wide format
#' set.seed(1)
#' dataW <- GenData(n=104)$d
#'
#' ## estimate mixed model
#' analyzeData(dataW, getinfo = FALSE)
#' analyzeData(dataW, getinfo = TRUE)
#' analyzeData(dataW, ddf = "satterthwaite", getinfo = TRUE)
#'


## * analyzeData (code)
#' @export
analyzeData <- function(d, ddf = "nlme", getinfo = TRUE, trace = TRUE){

    requireNamespace("nlme")

    ## ** normalize arguments
    ddf <- match.arg(ddf, choices = c("nlme","satterthwaite"))
    long <- wide2long(d, rm.na = TRUE) ## remove missing value but not pipeline value from the data.set
    ## head(long[order(long$id),],n=10)
    ## summary(long)

    ## ** fit gls model
    ctrl <- nlme::glsControl(opt='optim')

    m <- nlme::gls(X ~ baseline*visit + Z*visit,
                   data = long,
                   correlation = nlme::corSymm(form=~visit.num|id), 
                   weights = nlme::varIdent(form=~1|visit),
                   method = "REML",
                   na.action = stats::na.exclude)

    ## ** store estimate from the gls model    
    out <- list(d.long=long,
                fit=m,
                Info = 1/m$varBeta["Z1","Z1"], 
                estimate = as.double(stats::coef(m)["Z1"]),
                se = as.double(sqrt(m$varBeta["Z1","Z1"])),
                statistic = NA,
                df = NA,
                p.value = NA
                )
    out$statistic <- abs(out$estimate)/out$se

    ## computation of the p-value  
    if(ddf=="nlme"){ ## with weird estimator of the degree of freedom from nlme::gls 
        out$df <- as.double(m$dims$N - m$dims$p)
        out$p.value <- as.double(2*(1-stats::pt(out$statistic, df = out$df)))
        ## summary(m)$tTable["Z1","p-value"]
    }else{ ## or using satterthwaite approximation
        requireNamespace("emmeans")
        if(trace){
            groupTest <- emmeans::emmeans(m, specs = ~Z|visit, data = long[!is.na(long$X),,drop=FALSE])
        }else{
            groupTest <- suppressMessages(emmeans::emmeans(m, specs = ~Z|visit, data = long[!is.na(long$X),,drop=FALSE]))
        }
        e.satterthwaite <- summary(graphics::pairs(groupTest, reverse = TRUE), by = NULL, infer = TRUE, adjust = "none")
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
        if(any(abs(out$interim$vcov - stats::vcov(m))>1e-10)){
            warning("Something went wrong when extracting the variance covariance matrix from gls \n",
                    "largest discrepancy between getInformation and gls: ",max(abs(out$interim$vcov - stats::vcov(m))))
     
        }
    }

    ## ** Exportpp
    class(out) <- "lmmGSD"
    return(out)
}


## * wide2long
wide2long <- function(d, rm.na = FALSE, id.na = NULL, Z.id.na = TRUE){

    dW <- as.data.frame(d)

    ## get the number of timepoints
    col.time <- grep("^t([[:digit:]]+)$",names(dW), value = TRUE)
    col.X <- grep("^X([[:digit:]]+)$*",names(dW), value = TRUE)
    col.missing <- grep("^missing([[:digit:]]+)$*",names(dW), value = TRUE)
    if(length(col.missing)==0){
        dW$missing1 <- as.numeric(is.na(dW$X1))
        dW$missing2 <- as.numeric(is.na(dW$X2))
        dW$missing3 <- as.numeric(is.na(dW$X3))
        col.missing <- grep("^missing([[:digit:]]+)$*",names(dW), value = TRUE)
    }
    n.times <- length(grep("t",names(dW)))

    ## convert to long format
    dL <- stats::reshape(dW,
                         idvar = c("id","Z","missing1","t1","X1"),
                         timevar = "visit",
                         v.names = c("missing","time","X"),
                         varying = list(c("missing2","missing3"),
                                        c("t2","t3"),
                                        c("X2","X3")),
                         direction = "long")
    rownames(dL) <- NULL

    names(dL)[names(dL)=="X1"] <- "baseline"
    dL$id <- factor(dL$id)
    dL$Z <- factor(dL$Z)
    dL$visit.num <- as.numeric(dL$visit)
    dL$visit <- stats::relevel(as.factor(dL$visit),ref=n.times-1)## relevel to make the treatment coefficient have the interpretation of threatment effect at end of follow-up
    ## add missing values for not yet observed ids
    if(!is.null(id.na)){
        index.notObs <- which(dL$id %in% id.na)
        dL[index.notObs, c("missing1","missing")] <- -1
        dL[index.notObs, c("t1","baseline","time","X")] <- NA
        if(Z.id.na){
            dL[index.notObs, c("Z")] <- NA
        }
    }
    dL.save <- dL

    ## remove lines corresponding to missing data
    if(rm.na){
        dL <- dL[(dL$missing1<=0),] ## at baseline 
        dL <- dL[(dL$missing<=0),] ## or at follow-up
        dL$id <- droplevels(dL$id)
        dL$Z <- droplevels(dL$Z)
        dL$visit <- droplevels(dL$visit)
    }


    ## export
    attr(dL,"df.allobs") <- dL.save
    return(dL)
}
