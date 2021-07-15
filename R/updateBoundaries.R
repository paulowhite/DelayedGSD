## * updateBoundaries (documentation)
#' @title Update Boundaries of a GSD
#' @description Recompute the boundaries based according to the current information.
#'
#' @param x Object of type \code{delayedGSD}, typically output from \code{\link{CalcBoundaries}}.
#' @param lmm Optional argument used to update the information. Object of type \code{lmmGSD}, typically output from \code{\link{analyzeData}}.
#' @param k [integer] Index of the analysis.
#' @param type.k [character] Type of analysis: \code{"interim"} (after continuing recruitment),
#' \code{"decision"} (after stopping recruitment for efficacy or futility),
#' or \code{"final"} (after reaching the last stage of the trial).
#' @param trace [logical] should the execution of the function be traced?
#' @param update.stage [logical] should the arguments \code{k} and \code{type.k} be used to update to stage of the trial?
#' 

## * updateBoundaries (example)
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
#' ## generate data with all data for in case trial completes
#' set.seed(10)
#' theData <- GenData(n=theN*2,delta=theDelta*0.8,ar=5)  
#' 
#' theAR <- 10  #accrual rate (pt per month)
#' theDelay <- 0.7500001  #time in months to process data
#' tau.i <- theData$d$t3[theN + ceiling(theAR*theDelay)] #time point at which to do IA
#'
#'
#' #### Analyse data at the first interim ####
#' theInterimData <- SelectData(theData$d, t = tau.i, Delta.t = theDelay)
#' 
#' myLMM <- analyzeData(theInterimData)
#' myBound1 <- updateBoundaries(myBound0, lmm = myLMM, k = 1, type.k = "interim")
#' print(myBound1)
#' print(myBound1, planned = FALSE)
#' print(myBound1, planned = "only")
#'
#' par(mfrow = c(1,2))
#' plot(myBound1, planned = "only")
#' plot(myBound1)
#' 
#' #### Analyse data at the final ####
#' theFinalData <- SelectData(theData$d, t = 1e7, Delta.t = theDelay) 
#' 
#' myLMM <- analyzeData(theFinalData)
#' myBound2 <- updateBoundaries(myBound1, lmm = myLMM, k = 2, type.k = "final")
#' plot(myBound2)

## * updateBoundaries (code)
#' @export
updateBoundaries <- function(object, lmm = NULL, k, type.k, update.stage = TRUE, trace = TRUE){

    kMax <- object$kMax
    Info.max <- object$Info.max
    uk <- object$uk
    lk <- object$lk
    ck <- object$ck
    bindingFutility <- object$bindingFutility
    
    ## ** [optional] update information
    if(!is.null(lmm)){
        if(type.k %in% c("interim","final")){
            object$lmm[[k]] <- lmm
        }else{
            object$lmm[[k+1]] <- lmm
        }
        ## object$Info.i [1] 2.193777 3.640132
        object <- updateInformation(object, lmm = lmm, k = k, type.k = type.k, update.stage = FALSE)
    }

    Info.i <- object$Info.i
    Info.d <- object$Info.d
    
    ## ** update boundaries
    ## *** at interim
    if(type.k == "interim"){
     
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
                                        trace=trace,
                                        bindingFutility = bindingFutility)

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
    if(type.k == "decision"){

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
                                        trace=trace,
                                        bindingFutility = bindingFutility)

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
                                  bindingFutility=bindingFutility)
            object$ck[k:(kMax-1)]  <- c(newBound2[k], rep(NA,kMax-k-1))
        }else{
            ## What to do?
            stop("Do not know how to deal when Information decrease between interim and decision. \n")
        }
    }
    
    ## *** at final
    if(type.k == "final"){
        ## Should is be Info.max? or min(Info.i[k],Info.max)?
        ## Otherwise no need to update boundaries ...
        timing <- c(Info.i[1:(k-1)],object$Info.max)/object$Info.max

        if(bindingFutility){
            StandardDesign <- gsDesign(k=k, test.type=3, alpha=object$alpha, beta=object$beta,
                                       timing=timing,   #Need to double check that this is OK
                                       sfu=gsDesign::sfPower, sfupar=object$gammaA,
                                       sfl=gsDesign::sfPower, sflpar=object$gammaB)
            ## StandardDesign <- gsDesign(k=k, test.type=3, alpha=object$alpha, beta=object$beta,
            ##                            timing=timing,   #Need to double check that this is OK
            ##                            n.fix=1,n.I=Info.i,maxn.IPlan=object$InflationFactor,
            ##                            sfu=gsDesign::sfPower, sfupar=object$gammaA,
            ##                            sfl=gsDesign::sfPower, sflpar=object$gammaB)
        } else {
            StandardDesign <- gsDesign(k=k, test.type=4, alpha=object$alpha, beta=object$beta,
                                       timing=timing,   #Need to double check that this is OK
                                       n.I=Info.i, maxn.IPlan=object$InflationFactor,
                                       sfu=gsDesign::sfPower, sfupar=object$gammaA,
                                       sfl=gsDesign::sfPower, sflpar=object$gammaB)
            ## StandardDesign <- gsDesign(k=k, test.type=4,alpha=object$alpha,beta=object$beta,
            ##                            timing=timing,
            ##                            n.fix=1,n.I=Info.i,maxn.IPlan=onbject$InflationFactor,
            ##                            sfu=gsDesign::sfPower, sfupar=object$gammaA,
            ##                            sfl=gsDesign::sfPower, sflpar=object$gammaB)
        }

        object$lk[k]  <- StandardDesign$upper$bound[k]
        object$uk[k]  <- StandardDesign$upper$bound[k]
    }


    ## ** export stage
    if(update.stage){
        object$stage$k <- k
        object$stage$type <- type.k
    }

    ## ** export object
    return(object)
}
