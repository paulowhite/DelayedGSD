#' @examples
#' 
#' #### Planning #####
#' theAlpha <- 0.025
#' theBeta <- 0.2
#' theDelta <- 1.5
#' theK <- 2
#' theN <- 82
#' 
#' myBound <- CalcBoundaries(kMax=theK,
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
#' myBound <- update(myBound, data = theObsData, k = 1, analysis = "interim")
#'

## * update.delayedGSD (code)
#' @export
update.delayedGSD <- function(object, lmm, data, k, analysis, ...){

    ## ** check input arguments
    analysis <- match.arg(analysis, c("interim","decision"))

    ## check that we are at the next stage
    if(object$stage["k"]==0){
        if(k!=1){
            stop("Argument \'k\' should be 1. \n")
        }
        if(analysis!="interim"){
            stop("Argument \'k\' should be 1. \n")
        }
    }else{
    }
    object$stage["k"]
    
    ## ** get information
    ## fit mixed model and extract information
    if(missing(lmm) || is.null(lmm)){
        lmm <- AnalyzeData(data, getinfo = TRUE)
    }
    browser()
    lmm$getInformation["decision"]
    Ik=c(IA$Info,b1$Imax)
    
    ##     
    Bounds <- CalcBoundaries(kMax=object$kMax,
                             sided=object$sided,
                             alpha=object$alpha, 
                             beta=object$beta,  
                             InfoR.i=InfoR.i,  
                             gammaA=object$gammaA,
                             gammaB=object$gammaB,
                             method=object$method,
                             cNotBelowFixedc=object$cNotBelowFixedc, 
                             delta=object$delta,  
                             InfoR.d=InfoR.d,  
                             trace=FALSE,
                             bindingFutility = object$bindingFutility)

}
