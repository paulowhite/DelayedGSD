### getInformation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 11 2020 (10:18) 
## Version: 
## Last-Updated: dec 11 2020 (23:09) 
##           By: Brice Ozenne
##     Update #: 175
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * getInformation (documentation)
#' @title Extract information relative to a parameter.
#' @description Extract information relative to a parameter.
#'
#' @param object a \code{ttest} object or a \code{gls} object.
#' @param name.coef A character indicating relative to which coefficient the information should be computed.
#' @param method A character (\code{"direct"}, \code{"explicit"}, \code{"inflation"})
#' or a list containing the theoretical value of the variance-covariance matrix for each treatment group.
#' @param ... not used. For compatibility with the generic method.
#' 
#' @details Argument \bold{method}: \cr
#' \itemize{
#' \item \code{"direct"} use the square of the inverse of the estimated standard error for the regression coefficient to estimate the information.
#' \item \code{"explicit"} plug in the formula for the information / variance, the quantities estimated by the model (e.g. residual variance, effective sample size).
#' \item \code{"inflation"} same as \code{"explicit"} except it uses the number of samples instead of the number of samples with observations as the sample size when doing the plug-in.
#' \item \code{"pooling"} pool the current information with the future information of observations with missing data.
#' }
#' 
#' @export
`getInformation` <- function(object, method, ...) UseMethod("getInformation")

## * getInformation (examples)
#' @examples
#' library(nlme)
#' n <- 1e2
#'
#' #### Single endpoint ####
#' ## simulate data
#' set.seed(10)
#' X <- rnorm(n)
#' Y <- rnorm(n)
#' df <- rbind(data.frame(group=0,value=X),
#'             data.frame(group=1,value=Y))
#'
#' ## t-test
#' getInformation(ttest(value~1, df[df$group==0,])) ## only work for R>=4.0
#' getInformation(ttest(x = X)) ##  based on the estimate vcov
#' getInformation(ttest(X), variance = NULL) ## based on the estimated variance
#' getInformation(ttest(X), variance = 1) ## based on a variance of 1
#' getInformation(ttest(X), variance = 2) ## based on a variance of 1
#' getInformation(ttest(rnorm(length(X))), variance = 2) ## note: the X values do not matter here
#' 
#' getInformation(ttest(value~group, data = df), method = "direct")
#' getInformation(ttest(x = X, Y), variance = NULL)
#' getInformation(ttest(X,Y), variance = 1:2) ## information with a variance of 1 in one group and 2 in the other group
#' getInformation(ttest(X,Y), variance = c(1,3)) ## information with a variance of 1 in one group and 3 in the other group
#'
#' ## gls
#' library(nlme)
#'
#' e.gls <- gls(value~1, df[df$group==0,], na.action = na.omit)
#' getInformation(e.gls, name.coef = "(Intercept)", type = "prediction-pooling")
#' getInformation(e.gls, name.coef = "(Intercept)", method = "explicit")
#' getInformation(e.gls, name.coef = "(Intercept)", method = list(matrix(1,1,1)))
#' 
#' e.gls <- gls(value~group, data = df, weights = varIdent(form=~1|group))
#' getInformation(e.gls, name.coef = "(Intercept)", method = "direct")
#' getInformation(e.gls, name.coef = "(Intercept)", method = "explicit")
#' getInformation(e.gls, name.coef = "group", method = "direct")
#' getInformation(e.gls, name.coef = "group", method = "explicit")
#' getInformation(e.gls, name.coef = "group", method = list(matrix(1,1,1),matrix(2,1,1)))
#' 
#' #### Two endpoints ####
#' ## simulate data
#' library(mvtnorm)
#' 
#' set.seed(10)
#' X <- rmvnorm(n, sigma = 0.5 + diag(0.5,2))
#' Y <- rmvnorm(n, sigma = 0.5 + diag(0.5,2))
#' df <- rbind(data.frame(id = paste0("id",1:n),group=0,time=0,value=X[,1]),
#'             data.frame(id = paste0("id",1:n),group=0,time=1,value=X[,2]),
#'             data.frame(id = paste0("id",n+1:n),group=1,time=0,value=Y[,1]),
#'             data.frame(id = paste0("id",n+1:n),group=1,time=1,value=Y[,2]))
#' df <- df[order(df$id),]
#' 
#'
#' ## gls
#' e.gls <- gls(value~time-1, data = df[df$group==0,],
#'              correlation = corSymm(form=~1|id),
#'              weights = varIdent(form=~1|time))
#' getInformation(e.gls, name.coef = "time", method = "direct")
#' getInformation(e.gls, name.coef = "time", method = "explicit")
#' getInformation(e.gls, name.coef = "time", method = "pooling")
#' getInformation(e.gls, name.coef = "time", method = list(diag(1:2)))
#' 
#' e.gls <- gls(value~time*group, data = df,
#'              correlation = corSymm(form=~1|id),
#'              weights = varIdent(form=~1|time*group))
#' getInformation(e.gls, name.coef = "time:group", method = "direct")
#' getInformation(e.gls, name.coef = "time:group", method = "explicit")
#' getInformation(e.gls, name.coef = "time:group", method = "pooling")
#' ls.Sigma <- list(getVarCov(e.gls, individual = "id1"), getVarCov(e.gls, individual = "id101"))
#' getInformation(e.gls, name.coef = "time:group", method = ls.Sigma)
#' getInformation(e.gls, name.coef = "time:group", method = "inflation")
#' 

## * getInformation.matrix
getInformation.matrix <- function(object, variance, variance.group, ...){

    if(!is.list(variance)){variance <- list(variance)}
    
    browser()
    for(iId in 1:n.id){ ## iId <- 1
        if(is.null(object$modelStruct$corStruct) && is.null(object$modelStruct$varStruct)){
            iX <- X[iId,,drop=FALSE]
            iSigma <- method[[1]]
        }else if(!is.null(object$modelStruct$corStruct)){
            iIndex <- which(object$group==vec.id[iId])
            iX <- X[iIndex,,drop=FALSE]
            iSigma <- method[[unique(index.Sigma[iIndex])]]
        }else if(!is.null(object$modelStruct$varStruct)){
            iX <- X[iId,,drop=FALSE]
            iSigma <- method[[index.Sigma[iId]]]
        }
        Info <- Info + t(iX) %*% solve(iSigma) %*% iX
    }
        
    return(1/solve(Info)[name.coef,name.coef])
}

## * getInformation.htest
getInformation.ttest <- function(object, type = "estimation", variance, ...){

    ## ** normalize arguments
    type <- match.arg(type, c("estimation","prediction"))
    match.arg(object$method, c("One Sample t-test","Welch Two Sample t-test"))

    ## ** compute standard error for the mean/difference in mean
    if(object$method=="One Sample t-test"){
        z <- list(object$args$x)
    }else if(object$method=="Welch Two Sample t-test"){
        z <- list(x=object$args$x, y=object$args$y)
    }

    if(missing(variance)){
        
        if(type=="estimation"){ ## output the current estimated information
        
            se <- object$stderr
        
        }else if(type=="prediction"){ ## output the predicted information if there was no missing value

            se <- sqrt(Reduce("+",lapply(z, function(iZ){var(iZ, na.rm = TRUE)/length(iZ)})))

        }
        
    }else{

        ## get n
        if(type=="estimation"){
            n <- sapply(z, function(iZ){sum(!is.na(iZ))})
        }else if(type == "prediction"){
            n <- sapply(z, function(iZ){length(iZ)})
        }
        
        ## get variance
        if(!is.null(variance)){
            if(is.list(variance)){variance <- unlist(variance)}
            if(object$method=="One Sample t-test"){
                if(length(variance) != 1){
                    stop("The number of variance parameters contained in the \'method\' argument do not match the number of groups.\n")
                }
            }else if(object$method=="Welch Two Sample t-test"){
                if(length(variance) != 2){
                    stop("The number of variance parameters contained in the \'method\' argument do not match the number of groups.\n")
                }
            }
        }else{
            variance <- sapply(z, function(iZ){var(iZ, na.rm = TRUE)})
        }
        
        se <- sqrt(sum(variance/n))
        
    }
        
    ## ** export information
    return(as.double(1/se^2))
}

## * getInformation.gls
getInformation.gls <- function(object, name.coef, type = "estimation", variance, data = NULL,...){

    ## ** normalize arguments
    type <- match.arg(type, c("estimation","prediction-inflation","prediction-pooling"))

    ## if(!is.null(object$na.action)){
    ##     stop("Argument \'na.action\' should be null when calling gls. \n",
    ##          "The rows with NA should be excluded in the dataset by the user not by gls. \n")
    ## }
    if(is.null(data)){
        data <- try(nlme::getData(object), silent = TRUE)
        if(inherits(data,"try-error")){
            stop("Could not retrieve the data used to fit the gls model. \n",
                 "Consider passing the data via the \'data\' argument. \n")
        }
    }

    name.allcoef <- names(coef(object))
    n.allcoef <- length(name.allcoef)
    name.coef <- match.arg(name.coef, name.allcoef)

    ## ** compute standard error for the mean/difference in mean
    if(missing(variance)){

        if(type == "estimation"){
            
            se <- 1/vcov(object)[name.coef,name.coef]
            
        }else if(type== "prediction-pooling"){

            ## extract current variance estimate
            var.lmm <- vcov(object)[name.coef,name.coef]

            ## find number of observation per cluster
            if(is.null(object$modelStruct$corStruct)){
                n.interim.full <- object$dims$N
                n.interim.proxy <- length(object$na.action)
                n.decision <- n.interim.full + n.interim.proxy

                rho <- 0
                var.cc <- var.lmm
            }else{
                n.obs.group <- sapply(attr(object$modelStruct$corStruct,"covariate"), length)
                n.decision <- length(unique(object$group))
                n.interim.proxy <- sum(n.obs.group==1)
                n.interim.full <- sum(n.obs.group==2)
                rho <- coef(object$modelStruct$corStruct, unconstrained = FALSE)
                if(any(n.obs.group %in% 1:2 == FALSE)){
                    stop("Can only deal datasets with 1 or 2 observations per cluster \n")
                }
                object$groups %in% names(which(n.obs.group==2))

                ## complete case analysis (should be equivalent to the t-test)
                object.cc <- update(object, data = data[index.cc,,drop=FALSE])
                var.cc <- vcov(object.cc)[name.coef,name.coef]
            }
            
            ## se
            se <- sqrt(var.lmm - (1-rho^2) * var.cc * (n.interim.proxy/n.decision))
        
    }else if(method == "inflation"){
        if(is.null(object$modelStruct$corStruct)){
            stop("Inflation only makes sense in presence of correlated endpoints. \n",
                 "No correlation structure could be found in the object. \n")
        }

        n.obs.group <- sapply(attr(object$modelStruct$corStruct,"covariate"), length)
        if(any(n.obs.group %in% 1:2 == FALSE)){
            stop("Can only deal datasets with 1 or 2 observations per cluster \n")
        }

        var.lmm <- vcov(object)[name.coef,name.coef]
        rho <- coef(object$modelStruct$corStruct, unconstrained = FALSE)
        n.decision <- length(unique(object$group))
        n.interim.proxy <- sum(n.obs.group==1)
        n.interim.full <- sum(n.obs.group==2)
        
        ## complete case analysis (should be equivalent to the t-test)
        data.cc <- data[object$groups %in% names(which(n.obs.group==2)),]
        object.cc <- update(object, data = data.cc)
        vec.id2 <- object.cc$groups
        
        if(is.null(object.cc$modelStruct$varStruct)){
            if(is.null(vec.id2)){vec.id2 <- 1:NROW(data.cc)}
            ls.sigma <- list(unclass(getVarCov(object.cc, individual = vec.id2[1])))
        }else{
            Sigma.pattern <- unique(do.call(rbind,tapply(attr(object.cc$modelStruct$varStruct,"groups"),object.cc$groups, function(iVec){list(iVec)})))
            ls.sigma <- lapply(rownames(Sigma.pattern), function(iID){unclass(getVarCov(object.cc, individual = iID))})
        }
        info.cc <- getInformation(e.gls, name.coef = name.coef, method = ls.sigma)
        return(info.cc*n.decision/n.interim.full)
        
    }
    }else{
    }
    ## else{

    ##     ## identify independent clusters of observations
    ##     vec.id <- unique(object$group)
    ##     if(length(vec.id)==0){
    ##         vec.id <- 1:NROW(data)
    ##     }
    ## }
    
    

    ## ## ** prepare
    ## sigma2 <- stats::sigma(object)^2
    ## if(!is.null(object$modelStruct$varStruct)){
    ##     vec.sigma2 <- setNames(c(1, coef(object$modelStruct$varStruct, unconstrained = FALSE)^2) * sigma2,attr(object$modelStruct$varStruct,"groupNames"))
    ## }

    ## ## ** compute information
    ## if(is.list(method)){
    ##     X <- model.matrix(formula(object), data = data)
    ##     if(is.null(object$modelStruct$corStruct)){
    ##         n.id <- NROW(X)
    ##     }else{
    ##         n.id <- length(vec.id)
    ##     }
        
    ##     Info <- matrix(0, nrow = n.allcoef, ncol = n.allcoef,
    ##                    dimnames = list(name.allcoef, name.allcoef))
 
    ##     if(is.null(object$modelStruct$varStruct)){
    ##         index.Sigma <- rep(1,length(vec.id))
    ##     }else if(!is.null(object$modelStruct$corStruct)){
    ##         Sigma.pattern <- unique(do.call(rbind,tapply(attr(object$modelStruct$varStruct,"groups"),object$groups, function(iVec){list(iVec)})))
    ##         index.Sigma <- apply(apply(Sigma.pattern, 1, function(iPattern){attr(object$modelStruct$varStruct,"groups") %in% iPattern}),1,which)
    ##     }else{
    ##         index.Sigma <- as.numeric(as.factor(attr(object$modelStruct$varStruct,"groups")))
    ##     }
    ##     if(length(unique(sort(index.Sigma))) != length(method) || any(unique(sort(index.Sigma)) %in% 1:length(method) == FALSE)){
    ##         stop("The number of covariance matrices contained in the \'method\' argument do not match the number of covariate levels.\n")
    ##     }

    ##     getInformation(X, method = method, variance.group  )
    ##     for(iId in 1:n.id){ ## iId <- 1
    ##         if(is.null(object$modelStruct$corStruct) && is.null(object$modelStruct$varStruct)){
    ##             iX <- X[iId,,drop=FALSE]
    ##             iSigma <- method[[1]]
    ##         }else if(!is.null(object$modelStruct$corStruct)){
    ##             iIndex <- which(object$group==vec.id[iId])
    ##             iX <- X[iIndex,,drop=FALSE]
    ##             iSigma <- method[[unique(index.Sigma[iIndex])]]
    ##         }else if(!is.null(object$modelStruct$varStruct)){
    ##             iX <- X[iId,,drop=FALSE]
    ##             iSigma <- method[[index.Sigma[iId]]]
    ##         }
    ##         Info <- Info + t(iX) %*% solve(iSigma) %*% iX
    ##     }
        
    ##     return(1/solve(Info)[name.coef,name.coef])

    ## }else if(method=="explicit"){

    ##     X <- model.matrix(formula(object), data = data)
    ##     if(is.null(object$modelStruct$corStruct)){
    ##         n.id <- NROW(X)
    ##     }else{
    ##         n.id <- length(vec.id)
    ##     }
        
    ##     Info <- matrix(0, nrow = n.allcoef, ncol = n.allcoef,
    ##                    dimnames = list(name.allcoef, name.allcoef))

    ##     if(any(diff(as.numeric(object$groups)) %in% c(0,1) == FALSE)){
    ##         stop("Dataset must be ordered by id for getVarCov to work properly \n")
    ##     }
        
    ##     for(iId in 1:n.id){ ## iId <- 1
            
    ##         if(is.null(object$modelStruct$corStruct)){
    ##             iX <- X[iId,,drop=FALSE]
    ##             if(is.null(object$modelStruct$varStruct)){
    ##                 iSigma <- matrix(sigma2)
    ##             }else{
    ##                 iSigma <- vec.sigma2[attr(object$modelStruct$varStruct,"groups")[iId]]
    ##             }
    ##         }else{
    ##             iIndex <- which(object$group==vec.id[iId])
    ##             iX <- X[iIndex,,drop=FALSE]
    ##             iSigma <- unclass(getVarCov(object, individual = vec.id[iId])) ## incorrect variance structure when using getVarCov
    ##         }
    ##         Info <- Info + t(iX) %*% solve(iSigma) %*% iX
    ##     }

    ##     return(1/solve(Info)[name.coef,name.coef])

    ## }
}


######################################################################
### getInformation.R ends here
