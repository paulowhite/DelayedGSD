#' @examples
#' 
#' #### Planning #####
#' theAlpha <- 0.025
#' theBeta <- 0.2
#' theDelta <- 1.5
#' theK <- 2
#' theN <- 82
#' 
#' myBound0 <- CalcBoundaries(kMax=theK,
#'                           sided=1,
#'                      alpha=theAlpha,
#'                      beta=theBeta,
#'                      InfoR.i=c(0.5,1),
#'                      gammaA=2,
#'                      gammaB=2,
#'                      method=1,
#'                      delta=theDelta,
#'                      InfoR.d=0.55)
#'
#' #### Simulate data ####
#' set.seed(10)
#' theData <- GenData(n=theN*2,delta=theDelta*0.8,ar=5)  #generate data with all data for in case trial completes
#' 
#' theAR <- 10  #accrual rate (pt per month)
#' theDelay <- 0.7500001  #time in months to process data
#' tau.i <- theData$d$t3[theN + ceiling(theAR*theDelay)] #time point at which to do IA
#'
#' theObsData <- SelectData(theData$d, t = tau.i, Delta.t = theDelay)  #data at IA when deciding whether to continue recruitment
#'
#' #### Analyse data at interim ####
#' myBound1 <- update(myBound0, data = theObsData, k = 1, analysis = "interim")
#' print(myBound1)
#' print(myBound1, planned = FALSE)
#' print(myBound1, planned = "only")
#'
#' par(mfrow = c(1,2))
#' plot(myBound1, planned = "only")
#' plot(myBound1)

## * update.delayedGSD (code)
#' @export
update.delayedGSD <- function(object, lmm, data, k, analysis, ...){

    ## ** check input arguments
    kMax <- object$kMax
    if(k == kMax){
        analysis <- "interim"
    }
    analysis <- match.arg(analysis, c("interim","decision"))

    ## check that we are at the next stage
    object.k <- object$stage$k
    object.decision <- object$stage$decision
    if(object.k==kMax){
        previous.stage <- "final"
    }else if(object.decision>0){
        previous.stage <- "decision"
    }else{
        previous.stage  <- "interim"
    }

    if(object.k==0){
        stage <- "interim"
        object.conclusionInterim <- as.character(NA)
        if(k!=1){
            stop("Argument \'k\' should be 1 just after the planning stage. \n")
        }
        if(analysis!="interim"){
            stop("Argument \'analysis\' should be \"interim\" just after the planning stage. \n")
        }
    }else{
        object.conclusionInterim <- object$conclusion["interim",object.k]
        if(is.na(object.conclusionInterim)){
            stop("Decision analysis must be made for stage ",object.k," before updating the boundaries. \n")
        }

        if(previous.stage %in% "decision"){
            stop("No more boundary to update when a decision analysis has been performed. \n")
        }
        if(previous.stage %in% "final"){
            stop("No more boundary to update when the final analysis has been performed. \n")
        }
        if(object.conclusionInterim=="continue"){
            if(k!=(object.k+1)){
                stop("Argument \'k\' should be ",object.k+1," after continuing recruitment following the interim of stage ",object.k,". \n")
            }
            if(k==kMax){
                stage <- "final"
            }else{
                stage <- "interim"
            }
            if(analysis!="interim"){
                stop("Argument \'analysis\' should be \"interim\" after continuing recruitment following the interim of stage ",object.k,". \n")
            }
        }else{
            stage <- "decision"
            if(k!=object.k){
                stop("Argument \'k\' should be ",object.k," after stopping recruitment following the interim of stage ",object.k,". \n")
            }
            if(analysis!="decision"){
                stop("Argument \'analysis\' should be \"decision\" after stopping recruitment following the interim of stage ",object.k,". \n")
            }
        }
    }
    
    ## ** update lmm
    ## fit mixed model and extract information
    if(missing(lmm) || is.null(lmm)){
        lmm <- AnalyzeData(data, getinfo = TRUE)
    }
    if(stage %in% c("interim","final")){
        object$lmm[[k]] <- lmm
    }else{
        object$lmm[[k+1]] <- lmm
    }

    ## ** update information
    if(stage %in% c("interim","final")){
        object$Info.i[k] <- as.double(lmm$getInformation["interim"])
    }
    if(stage %in% c("interim","decision")){
        ## update decision (even when doing interim) to ensure monotone information
        object$Info.d[k] <- as.double(lmm$getInformation["decision"])
    }

    Info.max <- object$Info.max
    Info.i <- object$Info.i
    Info.d <- object$Info.d
    uk <- object$uk
    lk <- object$lk
    ck <- object$ck
    
    ## ** update boundaries
    ## *** at interim
    if(stage == "interim"){
     
        if(Info.d[k] < Info.max){
            newBounds <- CalcBoundaries(kMax=kMax,
                                        sided=object$sided,
                                        alpha=object$alpha, 
                                        beta=object$beta,  
                                        InfoR.i=Info.i/object$Info.max,  
                                        gammaA=object$gammaA,
                                        gammaB=object$gammaB,
                                        method=object$method,
                                        cNotBelowFixedc=object$cNotBelowFixedc, 
                                        delta=object$delta,  
                                        InfoR.d=Info.d/object$Info.max,  
                                        trace=FALSE,
                                        bindingFutility = object$bindingFutility)

            object$lk  <- newBounds$lk
            object$uk  <- newBounds$uk
            object$ck  <- newBounds$ck
        }else{
            object$lk[k]  <- Inf
            object$uk[k]  <- -Inf
            object$ck[k]  <- NA
        }

    }

    ## *** at decision
    if(previous.stage == "decision"){

        if(Info.d[k] < Info.max){
            newBounds <- CalcBoundaries(kMax=kMax,
                                        sided=object$sided,
                                        alpha=object$alpha, 
                                        beta=object$beta,  
                                        InfoR.i=Info.i/object$Info.max,  
                                        gammaA=object$gammaA,
                                        gammaB=object$gammaB,
                                        method=object$method,
                                        cNotBelowFixedc=object$cNotBelowFixedc, 
                                        delta=object$delta,  
                                        InfoR.d=Info.d/object$Info.max,  
                                        trace=FALSE,
                                        bindingFutility = object$bindingFutility)

            object$lk  <- newBounds$lk
            object$uk  <- newBounds$uk
            object$ck  <- newBounds$ck

        }else if(object$conclusion["reason.interim",k]=="Imax reached"){
            ## recompute the decision boundary to spend all the alpha
            newBounds2 <- Method1(uk=uk[1:k],
                                  lk=lk[1:k],
                                  Info.i=Info.i[1:k],
                                  Info.d=Info.d[k],
                                  Info.max=Info.max,
                                  sided=object$sided,
                                  ImaxAnticipated=TRUE,
                                  rho=object$gammaA,
                                  alpha=object$alpha,
                                  bindingFutility=object$bindingFutility)
            object$ck[k:(kMax-1)]  <- c(newBound2[k], rep(NA,kMax-k-1))
        }else{
            ## What to do?
            stop("Do not know how to deal when Information decrease between interim and decision. \n")
        }
    }
    
    ## *** at final
    if(previous.stage == "final"){

        if(bindingFutility){
            StandardDesign <- gsDesign(k=k, test.type=3,alpha=object$alpha,beta=object$beta,
                                       timing=c(Info.i[1:(k-1)],object$Info.max)/object$Info.max,   #Need to double check that this is OK
                                       n.fix=1,n.I=Info.i,maxn.IPlan=object$InflationFactor,
                                       sfu=sfPower,sfupar=object$gammaA,
                                       sfl=sfPower,sflpar=object$gammaB)
        } else {
            StandardDesign <- gsDesign(k=k, test.type=4,alpha=object$alpha,beta=object$beta,
                                       timing=c(Info.i[1:(k-1)],object$Info.max)/object$Info.max,
                                       n.fix=1,n.I=Info.i,maxn.IPlan=object$InflationFactor,
                                       sfu=sfPower,sfupar=object$gammaA,
                                       sfl=sfPower,sflpar=object$gammaB)
        }

        object$lk[k]  <- StandardDesign$upper$bound[k]
        object$uk[k]  <- StandardDesign$upper$bound[k]
    }


    ## ** export object
    object$stage["k"] <- k
    object$stage["decision"] <- (analysis=="decision")
    return(object)
}
