### getInformation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 11 2020 (10:18) 
## Version: 
## Last-Updated: feb  9 2021 (17:02) 
##           By: Brice Ozenne
##     Update #: 690
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
#' @param name.coef [character] For which coefficient the information should be computed (only relevant for gls models).
#' @param type [character] Should the information be estimated only for the observed data (\code{"estimation"}, excluding missing values)
#' or for the full dataset  (\code{"prediction"}, including missing values).
#' @param method.prediction [character] method used to estimate the information relative to the missing values (only relevant for gls models).
#' Can be by increasing the information proportionally to the amount of missing values (\code{"inflation"})
#' or by approximating the variance of the estimator by a weighted average between the complete case estimator and the full information estimator (\code{"pooling"}).
#' or a list containing the theoretical value of the variance-covariance matrix for each treatment group.
#' @param variance [NULL, list] If \code{NULL} then the estimated variance-covariance matrix from the object is used.
#' Otherwise it should be a list containing the theoretical value of the variance-covariance matrix for each treatment group.
#' @param data [data.frame] The dataset relative to which the information should be computed. See details section.
#' @param details [logical] Should intermediate results be output. See details section.
#' @param ... not used. For compatibility with the generic details.

#' @method Argument \bold{data}: the dataset may contain missing value in the outcome but no in the covariates. Missing value in the outcome indicates that the information is not available at the interim analysis but will be come available at decision.
#'
#'
#' Argument \bold{details}: when using gls models, an attribute detail is added to the output which contain a list:\itemize{
#' \item decision: information at decision analysis
#' \item interim: information at the interim analysis using all available data
#' \item interim.cc: information at the interim analysis using a complete case analysis
#' \item n: sample size at decision, interim with complete observation, interim with only partial observations
#' }
#' 
#' @export
`getInformation` <- function(object, ...) UseMethod("getInformation")

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
#' getInformation(ttest(x = X)) ##  based on the estimate variance
#' getInformation(ttest(X), variance = 1) ## based on a variance of 1
#' getInformation(ttest(X), variance = 2) ## based on a variance of 2
#' getInformation(ttest(rnorm(length(X))), variance = 2) ## note: the X values do not matter here
#' 
#' getInformation(ttest(value~group, data = df))
#' getInformation(ttest(x = X, Y))
#' getInformation(ttest(X,Y), variance = 1:2) ## information with a variance of 1 in one group and 2 in the other group
#' getInformation(ttest(X,Y), variance = c(1,3)) ## information with a variance of 1 in one group and 3 in the other group
#'
#' ## gls
#' library(nlme)
#'
#' e.gls <- gls(value~1, data = df[df$group==0,])
#' getInformation(e.gls, name.coef = "(Intercept)")
#' getInformation(e.gls, name.coef = "(Intercept)", variance = matrix(1,1,1))
#' 
#' e.gls <- gls(value~group, data = df, weights = varIdent(form=~1|group))
#' getInformation(e.gls, name.coef = "group")
#' getInformation(e.gls, name.coef = "group", variance = list(matrix(1,1,1),matrix(2,1,1)))
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
#' 
#'
#' ## gls
#' e.gls <- gls(value~time-1, data = df[df$group==0,],
#'              correlation = corSymm(form=~1|id),
#'              weights = varIdent(form=~1|time))
#' getInformation(e.gls, name.coef = "time")  ## 1/vcov(e.gls)
#' getInformation(e.gls, name.coef = "time", variance = list(diag(1:2)))
#' 
#' e.gls <- gls(value~time*group, data = df,
#'              correlation = corSymm(form=~1|id),
#'              weights = varIdent(form=~1|time*group))
#' getInformation(e.gls, name.coef = "time:group") ## 1/vcov(e.gls)
#' getInformation(e.gls, name.coef = "time", variance = list(diag(1:2),diag(1:2)))
#' 

## * getInformation.matrix
getInformation.matrix <- function(object, variance, ...){
    ## object is the design matrix

    ## ** get hidden arguments
    dots <- list(...)
    index.variance <- dots$index.variance
    index.cluster <- dots$index.cluster

    ## ** prepare
    n.cluster <- length(index.variance)
    n.allcoef <- NCOL(object)
    name.allcoef <- colnames(object)
    Info <- matrix(0, nrow = n.allcoef, ncol = n.allcoef,
                   dimnames = list(name.allcoef, name.allcoef))
    
    ## ** compute information
    for(iId in 1:n.cluster){ ## iId <- 200
        iX <- object[index.cluster==iId,,drop=FALSE]
        iSigma <- variance[[index.variance[iId]]]
        Info <- Info + t(iX) %*% solve(iSigma) %*% iX
    }

    ## ** export
    return(Info)
}

## * getInformation.ttest
getInformation.ttest <- function(object, type = "estimation", variance = NULL, ...){

    ## ** normalize arguments
    type <- match.arg(type, c("estimation","prediction"))
    match.arg(object$method, c("One Sample t-test","Welch Two Sample t-test"))

    ## ** compute standard error for the mean/difference in mean
    if(object$method=="One Sample t-test"){
        z <- list(object$args$x)
    }else if(object$method=="Welch Two Sample t-test"){
        z <- list(x=object$args$x, y=object$args$y)
    }

    if(is.null(variance)){
        
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
        
        se <- sqrt(sum(variance/n))        
    }
        
    ## ** export information
    return(as.double(1/se^2))
}

## * getInformation.gls
getInformation.gls <- function(object, name.coef, type = "estimation", method.prediction  = "inflation",
                               variance = NULL, data = NULL, details = FALSE,...){

    ## ** normalize arguments
    type <- match.arg(type, c("estimation","prediction"))
    method.prediction <- match.arg(method.prediction, c("inflation","pooling"))

    if(is.null(data)){
        if(type == "estimation"){
            data <- try(nlme::getData(object), silent = TRUE)
        }else{
            data <- try(eval(object$call$data), silent = TRUE)
        }
        if(inherits(data,"try-error")){
            stop("Could not retrieve the data used to fit the gls model. \n",
                 "Consider passing the data via the \'data\' argument. \n")
        }
    }
    data <- as.data.frame(data)
    
    if(!is.null(object$modelStruct$varStruct) && !inherits(object$modelStruct$varStruct, what = "varIdent")){
        stop("Only handles \"varIdent\" variance structures \n.")
    }
    if(!is.null(object$modelStruct$corStruct) && !inherits(object$modelStruct$corStruct, what = c("corSymm","corCompSymm"))){
        stop("Only handles \"corSymm\" and \"corCompSymm\" correlation structures \n.")
    }
    
    if(!is.null(variance)){
        if(!is.list(variance)){
            variance <- list(variance)
        }
        if(any(sapply(variance,is.matrix)==FALSE)){
            stop("Argument \"variance\" must contain a list of matrices. \n")
        }
        
    }

    ## ** extract elements from model
    ## *** formula
    f.gls <- stats::formula(object)
    cluster.var <- utils::tail(all.vars(stats::formula(object$modelStruct$corStruct)),1)
    
    ## *** coefficients
    name.allcoef <- names(coef(object))
    n.allcoef <- length(name.allcoef)
    name.coef <- match.arg(name.coef, name.allcoef)

    ## ** design matrix and covariance pattern
    ## *** at decision
    ## keep all observations despite missing values in the response
    X.decision <- model.matrix(f.gls, data = model.frame(formula = f.gls, data = data, na.action = na.pass))
    ## instead of 
    ## X <- model.matrix(stats::formula(object), data = data)
    ## to keep rows with missing data
    resPattern.decision <- .getPattern(object, data = data, variance = variance, name.coef = name.coef)

    if(any(is.na(X.decision))){
        stop("Cannot handle missing values in the design matrix. \n")
    }

    ## *** at interim: full information approach
    ## keep all observations without missing values in the response
    X.interim <- model.matrix(f.gls, data = data)

    index.interim <- which(!is.na(data[[all.vars(f.gls)[1]]]))
    if(any(abs(X.decision[index.interim,]-X.interim)>1e-10)){ ## Sanity check
        warning("Something went wrong when selecting the data at interim. \n")
    }
    data.interim <- data[index.interim,,drop=FALSE]
    resPattern.interim <- .getPattern(object, data = data.interim, variance = variance, name.coef = name.coef)

    ## number of observation per cluster at interim
    resPattern.interim$nobs.vargroup <- setNames(sapply(resPattern.interim$variance.vargroup,NCOL)[resPattern.interim$index.vargroup], names(resPattern.interim$index.vargroup)) 

    ## *** at interim: complete case approach
    if(resPattern.decision$rep.full==1){
        X.interim.cc <- X.interim
        resPattern.interim.cc <- resPattern.interim
    }else{
        ## keep all observations for which the cluster has no missing values in the response
        cluster.cc <- names(resPattern.interim$nobs.vargroup[resPattern.interim$nobs.vargroup==resPattern.decision$rep.full])
        index.cc <- which(data[[cluster.var]] %in% cluster.cc)
        data.interim.cc <- data[index.cc,,drop=FALSE]

        X.interim.cc <- model.matrix(f.gls, data = data.interim.cc)
        resPattern.interim.cc <- .getPattern(object, data = data.interim.cc, variance = variance, name.coef = name.coef)
    }    
    
    ## ** compute information
    ## *** at decision
    info.decision <- getInformation(X.decision,
                                    variance = resPattern.decision$variance.vargroup,
                                    index.variance = resPattern.decision$index.vargroup,
                                    index.cluster = resPattern.decision$index.cluster)
    var.decision <- solve(info.decision)[name.coef,name.coef]

    score.decision <- .getScore(X.decision,
                                variance = resPattern.decision$variance.vargroup,
                                residuals = data[[all.vars(f.gls)[1]]] - X.decision %*% coef(object),
                                index.variance = resPattern.decision$index.vargroup,
                                index.cluster = resPattern.decision$index.cluster)

    ## *** at interim: full information approach
    info.interim <- getInformation(X.interim,
                                   variance = resPattern.interim$variance.vargroup,
                                   index.variance = resPattern.interim$index.vargroup,
                                   index.cluster = resPattern.interim$index.cluster)
    var.interim <- solve(info.interim)[name.coef,name.coef]

    ## *** at interim: complete case approach
    info.interim.cc <- getInformation(X.interim.cc,
                                      variance = resPattern.interim.cc$variance.vargroup,
                                      index.variance = resPattern.interim.cc$index.vargroup,
                                      index.cluster = resPattern.interim.cc$index.cluster)
    var.interim.cc <- solve(info.interim.cc)[name.coef,name.coef]
    ## object.cc <- update(object, data = data.interim.cc)
    ## coef(object.cc)
    ## coef(object)
    
    ## ** sample size
    n.decision <- resPattern.decision$n.cluster ## number of patients at decision
    n.interim <- resPattern.interim$n.cluster ## number of patients at interim (with at least one observation proxy or outcome)
    n.interim.cc <- resPattern.interim.cc$n.cluster ## number of patients at interim with complete data

    ## ** predict information
    if(type == "estimation"){
        out <- 1/var.interim
    }else if(method.prediction == "inflation"){ 
        out <- 1/var.decision ## same as (1/var.interim.cc) * (n.decision/n.interim.cc)
    }else if(method.prediction == "pooling"){ 

        ## test number of timepoints
        if(length(resPattern.decision$rho)>1){
            stop("Argument \'method.prediction\' can only be \"pooling\" when there is only two timepoints. \n")
        }

        ## get the number of observations missing at interim for the parameter of interest but not missing for the proxy
        n.interim.proxy <- .nProxy(n.interim = n.interim, X.interim = X.interim, n.interim.cc = n.interim.cc,
                                   X.decision = X.decision, name.coef = name.coef, index.interim = index.interim,
                                   data = data, data.interim = data.interim, cluster.var = cluster.var)

        ##  pool variance and deduce information
        info.current <- 1/(var.interim - (1-resPattern.decision$rho^2) * var.interim.cc * n.interim.proxy/n.interim) ## information if we had the outcome (instead of proxy) for all interim patients
        out <- info.current * (n.decision/n.interim) ## information add the new patients (if any) at interim
    }

    ## ** export
    if(details){
        attr(out,"details") <- list(decision = list(information = 1/solve(info.decision)[name.coef,name.coef],
                                                    vcov = solve(info.decision),
                                                    pattern = resPattern.decision),
                                    interim = list(information = 1/solve(info.interim)[name.coef,name.coef],
                                                   vcov = solve(info.interim),
                                                   pattern = resPattern.interim),
                                    interim.cc = list(information = 1/solve(info.interim.cc)[name.coef,name.coef],
                                                      vcov = solve(info.interim.cc),
                                                      pattern = resPattern.interim),
                                    n = c(decision = n.decision, interim = n.interim, interim.cc = n.interim.cc)
                                    )
    }
    
    return(out)
}

## * .getPattern
.getPattern <- function(object, data, name.coef, variance = NULL){

    ## ** reinitialize correlation/variance structure according to data
    if(!is.null(object$modelStruct$corStruct)){
        corStructNew <- do.call(class(object$modelStruct$corStruct)[1],
                                args = list(form = stats::formula(object$modelStruct$corStruct),
                                            value = coef(object$modelStruct$corStruct, unconstrained = FALSE))
                                )
        corStructNew <- Initialize(corStructNew, data = data)
        corgroups <- getGroups(corStructNew)
    }
    if(!is.null(object$modelStruct$varStruct)){
        varStructNew <- do.call(class(object$modelStruct$varStruct)[1],
                                args = list(form = stats::formula(object$modelStruct$varStruct),
                                            value = coef(object$modelStruct$varStruct, unconstrained = FALSE))
                                )
        varStructNew <- Initialize(varStructNew, data = data)
        vargroups <- getGroups(varStructNew)
    }
    
    ## ** number and position of the clusters among the observations
    if(is.null(object$modelStruct$corStruct)){ ## no correlation structure
        index.cluster <- 1:NROW(data)
        n.cluster <- NROW(data)

    }else{ ## correlation structure
        index.cluster <- setNames(as.numeric(corgroups), corgroups)
        n.cluster <- attr(corStructNew,"Dim")$M

    }

    ## ** extract number and position of unique residual variance-covariance structures
    if(is.null(object$modelStruct$varStruct)){ ## no variance structure
        index.vargroup <- setNames(rep(1,n.cluster), sort(unique(names(index.cluster))))
        n.vargroup <-  1
        Sigma.pattern <- list("sigma") ## not needed

    }else if(is.null(object$modelStruct$corStruct)){ ## no correlation structure
        index.vargroup <- as.numeric(as.factor(vargroups))
        n.vargroup <-  max(index.vargroup)
        Sigma.pattern <- as.list(levels(as.factor(vargroups))) ## not needeed
        
    }else{ ## variance and correlation structure
        ## variance parameter within cluster
        variance.per.cluster <- tapply(X = vargroups, INDEX = corgroups, FUN = function(iVec){list(iVec)}) ## WARNING MAY MESS UP THE ORDER
        variance.per.cluster <- variance.per.cluster[levels(corgroups)] ## PROPERLY REORDER
            
        ## unique variance patterns
        Sigma.pattern <- unique(variance.per.cluster)
        names(Sigma.pattern) <- sapply(Sigma.pattern, paste, collapse = "|")
        n.vargroup <- length(Sigma.pattern)
    
        ## associate each cluster to a variance structure
        index.vargroup <- sapply(variance.per.cluster, function(x){
            as.double(which(unlist(lapply(Sigma.pattern, identical, x))))
        })

    }
    
    ## ** correlation coefficient
    if(!is.null(variance)){
        if(NCOL(variance[[1]])>1 && NROW(variance[[1]])>1){
            rho <- variance[[1]][lower.tri(variance[[1]])]
        }else{
            rho <- 0
        }
        
    }else if(is.null(object$modelStruct$corStruct)){
        rho <- 0
    }else{
        rho <- coef(object$modelStruct$corStruct, unconstrained = FALSE)
    }
    
    ## ** full residual variance-covariance matrix
    rep.full <- max(sapply(Sigma.pattern, length))
    Sigma.pattern.full <- Sigma.pattern[sapply(Sigma.pattern, length) == rep.full]
    names(Sigma.pattern.full) <- sapply(Sigma.pattern.full, paste, collapse = "|")

    variance <- .getResVcov(object, variance = variance, rho = rho, Sigma.pattern = Sigma.pattern.full, rep.full = rep.full)

    ## ** (missing data) residual variance-covariance matrix    
    if(is.null(object$modelStruct$varStruct) || is.null(object$modelStruct$corStruct)){ ## no variance or correlation structure
        variance.vargroup <- variance
    }else{ ## variance and correlation structure
        name.vargroup <- names(Sigma.pattern)
        variance.vargroup <- setNames(vector(mode = "list", length = n.vargroup), name.vargroup)
        MSigma.pattern.full <- do.call(rbind,Sigma.pattern.full)

        for(i.vargroup in 1:n.vargroup){ ## i.vargroup <- 2
            iname.vargroup <- name.vargroup[i.vargroup]
            if(iname.vargroup %in% names(Sigma.pattern.full)){
                variance.vargroup[[i.vargroup]] <- variance[[iname.vargroup]]
            }else{
                test <- t(apply(MSigma.pattern.full, 1, `%in%`, Sigma.pattern[[iname.vargroup]])) ## t() because apply(,1,) returns the transposed version
                index <- which.max(rowSums(test))
                variance.vargroup[[i.vargroup]] <- variance[[names(index)]][test[index,],test[index,],drop=FALSE]
            }
        }
    }

    ## ** export
    return(list(index.cluster = index.cluster, ## index of the cluster per observation
                index.vargroup = index.vargroup, ## index of the variance pattern per cluster
                n.cluster = n.cluster, ## number of clusters
                n.vargroup = n.vargroup, ## number of variance patterns
                variance.vargroup = variance.vargroup, ## variance-covariance matrix associated to each variance pattern
                rho =  as.double(rho),
                rep.full = rep.full)) ## number of observation per cluster when no missing data
}

## * .getResVcov
.getResVcov <- function(object, variance, rho, Sigma.pattern, rep.full){

    ## ** build or check variance covariance patterns for fully observed cluster
    if(is.null(variance)){
        if(is.null(object$modelStruct$varStruct) && is.null(object$modelStruct$corStruct)){
            variance <- list(matrix(sigma(object)^2, nrow = 1, ncol = 1))
        }else if(is.null(object$modelStruct$varStruct)){
            variance <- unclass(getVarCov(e.gls))
        }else if(is.null(object$modelStruct$corStruct)){
            variance <- as.list((sigma(object)*coef(object$modelStruct$varStruct, allCoef = TRUE, unconstrained = FALSE))^2)
        }else{
            vec.sigma <- sigma(object)*coef(object$modelStruct$varStruct, allCoef = TRUE, unconstrained = FALSE)
            Mcor <- matrix(1, nrow = rep.full, ncol = rep.full)
            Mcor[lower.tri(Mcor)] <- rho
            Mcor[upper.tri(Mcor)] <- t(Mcor)[upper.tri(Mcor)]
            variance <- setNames(lapply(Sigma.pattern, function(x){
                tcrossprod(vec.sigma[x]) * Mcor
            }), names(Sigma.pattern))
        }
    }else{

        ## check format and names
        if(length(variance) != length(Sigma.pattern)){
            stop("Incorrect \'variance\' argument (length = ",length(variance),"). \n",
                 length(Sigma.pattern)," variance-covariance patterns detected (levels: \"",paste0(names(Sigma.pattern), collapse="\" \""),"\"). \n")
        }
        if(!is.null(names(variance))){
            if(any(names(variance) %in% names(Sigma.pattern) == FALSE)){
                stop("Incorrect \'variance\' argument (names:  \"",paste0(names(variance), collapse="\" \""),"\"). \n",
                     "variance-covariance patterns detected (levels: \"",paste0(names(Sigma.pattern), collapse="\" \""),"\"). \n")
            }else{
                variance <- variance[names(Sigma.pattern)]
            }
        }else{
            names(variance) <- names(Sigma.pattern)
        }

        
        ## check correlation
        ls.rho <- lapply(variance[sapply(variance,length)>1], function(x){cov2cor(x)[lower.tri(x)]})
        if(length(ls.rho)>1){
            M.rho <- do.call(rbind, ls.rho)
            test <- apply(M.rho, MARGIN = 2, FUN = function(x){max(abs(diff(x)))})
            if(test>1e-10){
                stop("When specificy the argument \"variance\", the correlation should be the same for all matrices. \n")
            }
        }
    }

    return(variance)
}

## * .getScore
.getScore <- function(object, variance, residuals, index.variance, index.cluster){

    ## ** prepare
    n.obs <- length(index.cluster)
    n.cluster <- length(index.variance)
    n.allcoef <- NCOL(object)
    name.allcoef <- colnames(object)
    Score <- matrix(NA, nrow = n.obs, ncol = n.allcoef,
                    dimnames = list(NULL, name.allcoef))
    
    ## ** compute score
    for(iId in 1:n.cluster){ ## iId <- 7
        iResidual <- residuals[index.cluster==iId,,drop=FALSE]
        iX <- object[index.cluster==iId,,drop=FALSE]
        iSigma <- variance[[index.variance[iId]]]
        Score[index.cluster==iId,] <- sweep(solve(iSigma) %*% iX, FUN = "*", MARGIN = 1, STATS = iResidual)
    }

    ## ** export
    return(Score)
}

## * .nProxy
## get the number of observations missing at interim for the parameter of interest but not missing for the proxy
.nProxy <- function(n.interim, n.interim.cc,
                    X.decision, X.interim, name.coef, index.interim,
                    data, data.interim, cluster.var){

    ## ** if only one parameter, then any missing value is missing for the parameter of interest
    if(NCOL(X.decision)==1){
        return(n.interim - n.interim.cc)
    }
    ## else it could be that some "lines" of the design matrix are not useful to estimate the parameter of interest
    ## and are only useful to estimate nuisance parameters.
    
    ## ** do we have full information about the patients included at interim?
    ## if yes we are comparing at the cluster level so we can indeed substract the sample size
    index.interim.full <- which(data[[cluster.var]] %in% unique(data.interim[[cluster.var]]))
    if(length(setdiff(index.interim.full, index.interim))==0){
        return(n.interim - n.interim.cc)
    }
    ## else could be that the added lines when doing full information do not benefit the parameter of interest
    ## and therefore the sample size is the same as complete case analysis.
    
    ## ** Is having full data for each cluster at interim better than the observed data at interim for estimating the parameter of interest?
    ## i.e. does score of the parameter of interest contains a linear combination of the score of other parameters, so that their contribution to the score of the parameter of interest will always be 0.
    ## example Score = [Score_alpha & Score_beta] where Score_alpha = [S_1 \\ S_2] and Score_beta = [S_2 \\ 0] so sum(S_2) must be 0 and therefore the last observations do not contribute to Score_alpha
    X.new <- X.decision[setdiff(index.interim.full, index.interim),,drop=FALSE]
    scorefit <- lm.fit(x = X.new[,setdiff(colnames(X.new),name.coef),drop=FALSE], y = X.new[,name.coef])
    if(any(abs(scorefit$residuals)>1e-10)){ ## not a linear combination
        return(n.interim - n.interim.cc)
    }
    
    ## check whether,  at interim (full information), the score of the parameter of interest contains the score of the linear combination, via the design matrix
    combin <- na.omit(coef(scorefit))
    X.combin <- Reduce("+",lapply(1:length(combin), function(x){
        X.interim[,names(combin)[x],drop=FALSE]*combin[x]
    }))
    test.0 <- (X.combin-X.interim[, name.coef,drop=FALSE])[which(abs(X.combin)>1e-10)]

    if((all(abs(test.0)<1e-10))){
        return(0) ## number of patients without the outcome measurment and only the proxy measurement 
    }else{
        return(n.interim - n.interim.cc) ## number of patients without the outcome measurment and only the proxy measurement
    }
}

######################################################################
### getInformation.R ends here
