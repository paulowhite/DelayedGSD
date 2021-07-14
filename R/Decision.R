### Decision.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 28 2020 (14:40) 
## Version: 
## Last-Updated: Jul 14 2021 (10:56) 
##           By: Brice Ozenne
##     Update #: 93
#----------------------------------------------------------------------
## 
### Commentary: 
##  Fucntion that maps the statistical results at interim (or final) analysis into
##  decision to reject, or continue or stop inclusind subjects.
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Decision (documentation)
##' @title Evaluate decision at an interim, decision or final analysis of a group sequential design with delayed endpoints
##' @description Maps the statistical results at an interim, decision or final analysis into a decision regarding whether to stop recruitment (interim) or whether to reject the null hypothesis (decision/final). Stopping boundaries are updated based on observed information and correct p-values, confidence intervals and point estimates are given.
##' 
##' @param analysis_res object of class ?? from AnalyzeData
##' @param planned_bnds object of class ?? from CalcBoundaries
##' @param k the stage at which the decision is to be made
##' @param analysis is it an interim, decision or final analysis
##' @param Info.i the observed (where possible) or expected information at each interim and the final analysis
##' @param InfoR.d the expected or observed information ratio at each decision analysis
##' @param PositiveIsGood whether a positive effect is considered beneficial (TRUE/FALSE)
##' @param Trace whether to print some messages
##' @param bindingFutility whether to use binding futility boundaries (use TRUE for binding)
##' @param plot whether the updated boundaries and the result should be plotted
##' 
##' @details zz
##' @return ff
##' @author Paul Blanche
##' 
##' @examples
##'
##' #### Planning #####
##' theAlpha <- 0.025
##' theBeta <- 0.2
##' theDelta <- 1.5
##' theK <- 2
##' theN <- 82
##' 
##' b1 <- CalcBoundaries(kMax=theK,
##'                      sided=1,
##'                      alpha=theAlpha,
##'                      beta=theBeta,
##'                      InfoR.i=c(0.5,1),
##'                      gammaA=2,
##'                      gammaB=2,
##'                      method=1,
##'                      delta=theDelta,
##'                      InfoR.d=0.55)
##' plot(b1)
##'
##' #### Simulate data ####
##' set.seed(10)
##' theData <- GenData(n=theN*2,delta=theDelta*0.8,ar=5)  #generate data with all data for in case trial completes
##' 
##' theAR <- 10  #accrual rate (pt per month)
##' theDelay <- 0.7500001  #time in months to process data
##' tau.i <- theData$d$t3[theN + ceiling(theAR*theDelay)] #time point at which to do IA
##'
##' theObsData <- SelectData(theData$d, t = tau.i, Delta.t= theDelay)  #data at IA when deciding whether to continue recruitment
##'
##' #### Analyse data at interim ####
##' lmm.interim <- AnalyzeData(theObsData)
##' IA <- Decision(analysis_res = lmm.interim, planned_bnds = b1, k = 1, analysis = "interim") 



## * Decision (code)
##' @export
Decision <- function(analysis_res,  #results from AnalyzeData
                     planned_bnds,  #results from CalcBoundaries
                     k=1, #at which phase are we?
                     analysis="interim", #is it an interim or decision or final analysis
                     Info.i,  #all I_k (information) from first interim to final analysis (observed where possible)
                     InfoR.d=NULL,  #expected or observed information RATIO at each decision analysis
                     PositiveIsGood=TRUE, # whether positive effect is good (i.e. positive trial)
                     Trace=TRUE, # whether to print some messages
                     bindingFutility=TRUE,
                     plot=T){  #should the boundaries and results be plotted?
    
    ## ** check input arguments
    if(!inherits(analysis_res,"lmmGSD")){
        stop("Argument \'analysis_res\' should be a \"lmmGSD\" object. \n",
             "(typically output by AnalyzeData). \n")
    }
    if(!inherits(planned_bnds,"delayedGSD")){
        stop("Argument \'analysis_res\' should be a \"delayedGSD\" object. \n",
             "(typically output by CalcBoundaries). \n")
    }
    if(is.null(InfoR.d) & !((analysis=="interim" & planned_bnds$method==1) | analysis=="final")){
        stop("InfoR.d must be specified.")
    }
    if(is.null(InfoR.d) & analysis=="interim" & planned_bnds$method==1 & k==1){
        InfoR.d <- (Info.i[1]/planned_bnds$Info.max + 1)/2 #Rk: in this case it does not matter, this value is not used.
        ## print(paste0("InfoR.d=",InfoR.d))
    }

    #update boundaries
    if(analysis!="final"){
        Bounds <- CalcBoundaries(kMax=planned_bnds$kMax,
                                 sided=planned_bnds$sided,  #one or two-sided
                                 alpha=planned_bnds$alpha,  #type I error
                                 beta=planned_bnds$beta,  #type II error
                                 InfoR.i=Info.i/planned_bnds$Info.max,  #planned or observed information rates
                                 gammaA=planned_bnds$gammaA,  #rho parameter for alpha error spending function
                                 gammaB=planned_bnds$gammaB,  #rho parameter for beta error spending function
                                 method=planned_bnds$method,  #use method 1 or 2 from paper H&J
                                 cNotBelowFixedc=planned_bnds$cNotBelowFixedc, # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                                 delta=planned_bnds$delta,  #effect that the study is powered for
                                 InfoR.d=InfoR.d,  #(expected) information ratio at each decision analysis
                                 trace=FALSE,
                                 bindingFutility = bindingFutility)
    }else{  ####Is this really needed? Couldn't we just use the Bounds from above?
      if(bindingFutility){
        StandardDesign <- gsDesign(k=k, test.type=3,alpha=planned_bnds$alpha,beta=planned_bnds$beta,
                                   timing=Info.i/max(Info.i),   #Need to double check that this is OK
                                   sfu=sfPower,sfupar=planned_bnds$gammaA,
                                   sfl=sfPower,sflpar=planned_bnds$gammaB)
      } else {
        StandardDesign <- gsDesign(k=k, test.type=4,alpha=planned_bnds$alpha,beta=planned_bnds$beta,
                                   timing=Info.i/max(Info.i),
                                   sfu=sfPower,sfupar=planned_bnds$gammaA,
                                   sfl=sfPower,sflpar=planned_bnds$gammaB)
      }
    }

    #test statistic
    Z <- analysis_res$statistic
    
    if(!PositiveIsGood){
        Z <- -Z
        if(Trace){message("Negative effect is good: schange sign of Z")}
    }
    
    #evaluate decision
    if(analysis=="interim"){
        details <- c(u=Bounds$uk[k],l=Bounds$lk[k])
        if(Z>Bounds$uk[k]){
            decision="Efficacy"
        } else if (Z<Bounds$lk[k]){
            decision="Futility"
        } else {
            decision="Continue"
        }
    } else if(analysis=="decision"){
        ## print(paste0("c=",Bounds$ck[k]))
        details <- c(c=Bounds$ck[k])
        if(Z>Bounds$ck[k]){
            decision="Efficacy"
        } else {
            decision="Futility"
        }
    } else if(analysis=="final"){
        #details <- c(critical=StandardDesign$criticalValues[k])
        details <- c(critical=StandardDesign$upper$bound[k])
        #if(Z>StandardDesign$criticalValues[k]){
        if(Z>StandardDesign$upper$bound[k]){
            decision="Efficacy"
        } else {
            decision="Futility"
        }
    } else {
        stop("Please specify interim or decision of final for analysis type")
    }
    
    if(plot){
        if(analysis=="final"){
            warning("Cannot produce the plot yet, when analysis=final.")
        }else{
            PlotBoundaries(Bounds)
            points(analysis_res$Info/planned_bnds$Info.max,Z,pch=4,cex=1.5)  #at some point we may wish to add that all previous analyses are plotted as well
        }
    }
    
    out <- list(decision=decision,details=details)
    out
}

#would be nice if we can have a predict_info function to avoid that InfoR.d needs to be an argument


## * Decision2 (documentation)
Decision2 <- function(analysis_res,  #results from AnalyzeData
                     planned_bnds,  #results from CalcBoundaries
                     k=1, #at which phase are we?
                     analysis="interim", #is it an interim or decision or final analysis
                     Info.i,  #all I_k (information) from first interim to final analysis (observed where possible)
                     InfoR.d,  #expected or observed information RATIO at each decision analysis
                     PositiveIsGood=TRUE, # whether positive effect is good (i.e. positive trial)
                     Trace=TRUE, # whether to print some messages
                     bindingFutility=TRUE,  #whether to use binding futility rules
                     ImaxAnticipated=FALSE,   #whether interim analysis k was skipped due to anticipated Imax reached
                     plot=T){  #should the boundaries and results be plotted?
  
  #check input arguments
  #if(is.null(InfoR.d) & !((analysis=="interim" & planned_bnds$method==1) | analysis=="final")){
  #  stop("InfoR.d must be specified.")
  #}
  #if(is.null(InfoR.d) & analysis=="interim" & planned_bnds$method==1 & k==1){
  #  InfoR.d <- (Info.i[1]/planned_bnds$Info.max + 1)/2 #Rk: in this case it does not matter, this value is not used.
  #  ## print(paste0("InfoR.d=",InfoR.d))
  #}
  
  #test statistic
  Z <- analysis_res$statistic
  
  if(!PositiveIsGood){
    Z <- -Z
    if(Trace){message("Negative effect is good: schange sign of Z")}
  }
  
  #decision at interim
  if(analysis=="interim"){
    #skip interim analysis if information is decreasing
    if(Info.i[k] < Info.i[k-1]){   
      decision <- "Continue"
      reason <- "Decreasing information"
      details <- c(u=NA,l=NA)
    #skip interim analysis and continue straight to decision analysis if it is anticipated that Imax will be reached at decision analysis k
    } else if(Info.i[k] > planned_bnds$Info.max | InfoR.d[k] > 1){
      decision <- "Stop recruitment"
      reason <- "Imax reached at interim or decision"
      details <- c(u=NA,l=NA)
    #update boundaries and evaluate decision as usual otherwise
    } else {
      Bounds <- CalcBoundaries(kMax=planned_bnds$kMax,
                               sided=planned_bnds$sided,  #one or two-sided
                               alpha=planned_bnds$alpha,  #type I error
                               beta=planned_bnds$beta,  #type II error
                               InfoR.i=Info.i/planned_bnds$Info.max,  #planned or observed information rates
                               gammaA=planned_bnds$gammaA,  #rho parameter for alpha error spending function
                               gammaB=planned_bnds$gammaB,  #rho parameter for beta error spending function
                               method=planned_bnds$method,  #use method 1 or 2 from paper H&J
                               cNotBelowFixedc=planned_bnds$cNotBelowFixedc, # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                               delta=planned_bnds$delta,  #effect that the study is powered for
                               InfoR.d=InfoR.d,  #(expected) information ratio at each decision analysis
                               trace=FALSE,
                               bindingFutility = bindingFutility)
      
      details <- c(u=Bounds$uk[k],l=Bounds$lk[k])
      if(Z>Bounds$uk[k]){
        decision <- "Stop recruitment"
        reason <- "Efficacy"
      } else if (Z<Bounds$lk[k]){
        decision <- "Stop recruitment"
        reason <- "Futility"
      } else {
        decision <- "Continue"
        reason <- "No boundaries crossed"
      }
    }
    
  #decision at decision analysis
  } else if(analysis=="decision"){
    #if information has decreased since interim analysis
    if(InfoR.d[k]*planned_bnds$Info.max < Info.i[k]){
      stop("cannot handle case Info.d[k] < Info.i[k]")
    #if Imax has been reached or anticipated
    } else if(InfoR.d[k] >= 1 | ImaxAnticipated){
      if(InfoR.d[k] >= 1 & !ImaxAnticipated){
        stop("don't know what to do if Id > Imax but this was not anticipated")
      } else {
        Bounds <- CalcBoundaries(kMax=planned_bnds$kMax,
                                 sided=planned_bnds$sided,  #one or two-sided
                                 alpha=planned_bnds$alpha,  #type I error
                                 beta=planned_bnds$beta,  #type II error
                                 InfoR.i=Info.i/planned_bnds$Info.max,  #planned or observed information rates
                                 gammaA=planned_bnds$gammaA,  #rho parameter for alpha error spending function
                                 gammaB=planned_bnds$gammaB,  #rho parameter for beta error spending function
                                 method=planned_bnds$method,  #use method 1 or 2 from paper H&J
                                 cNotBelowFixedc=planned_bnds$cNotBelowFixedc, # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                                 delta=planned_bnds$delta,  #effect that the study is powered for
                                 InfoR.d=InfoR.d,  #(expected) information ratio at each decision analysis
                                 trace=FALSE,
                                 bindingFutility = bindingFutility)
        critval <- Method1(uk=Bounds$uk[1:k],
                           lk=Bounds$lk[1:k],
                           Info.i=Bounds$Info.i[1:k],
                           Info.d=Bounds$Info.d[k],
                           Info.max=Bounds$Info.max,
                           sided=Bounds$sided,
                           cMin=Bounds$cMin,
                           ImaxAnticipated=TRUE,
                           rho=Bounds$gammaA,
                           alpha=Bounds$alpha,
                           bindingFutility=Bounds$bindingFutility)
        
        details <- c(c=critval)
        if(Z>critval){
          decision="Efficacy"
          reason <- "Efficacy"
        } else {
          decision="Futility"
          reason <- "Futility"
        }
        
      }
    #in the usual case
    } else {
      Bounds <- CalcBoundaries(kMax=planned_bnds$kMax,
                               sided=planned_bnds$sided,  #one or two-sided
                               alpha=planned_bnds$alpha,  #type I error
                               beta=planned_bnds$beta,  #type II error
                               InfoR.i=Info.i/planned_bnds$Info.max,  #planned or observed information rates
                               gammaA=planned_bnds$gammaA,  #rho parameter for alpha error spending function
                               gammaB=planned_bnds$gammaB,  #rho parameter for beta error spending function
                               method=planned_bnds$method,  #use method 1 or 2 from paper H&J
                               cNotBelowFixedc=planned_bnds$cNotBelowFixedc, # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                               delta=planned_bnds$delta,  #effect that the study is powered for
                               InfoR.d=InfoR.d,  #(expected) information ratio at each decision analysis
                               trace=FALSE,
                               bindingFutility = bindingFutility)
      
      details <- c(c=Bounds$ck[k])
      if(Z>Bounds$ck[k]){
        decision="Efficacy"
        reason <- "Efficacy"
      } else {
        decision="Futility"
        reason <- "Futility"
      }
    }
    
  #decision at the final analysis  
  } else if(analysis=="final"){
    if(bindingFutility){
      StandardDesign <- gsDesign(k=k, test.type=3,alpha=planned_bnds$alpha,beta=planned_bnds$beta,
                                 timing=c(Info.i[1:(k-1)],planned_bnds$Info.max)/planned_bnds$Info.max,   #Need to double check that this is OK
                                 n.fix=1,n.I=Info.i,maxn.IPlan=planned_bnds$InflationFactor,
                                 sfu=sfPower,sfupar=planned_bnds$gammaA,
                                 sfl=sfPower,sflpar=planned_bnds$gammaB)
    } else {
      StandardDesign <- gsDesign(k=k, test.type=4,alpha=planned_bnds$alpha,beta=planned_bnds$beta,
                                 timing=c(Info.i[1:(k-1)],planned_bnds$Info.max)/planned_bnds$Info.max,
                                 n.fix=1,n.I=Info.i,maxn.IPlan=planned_bnds$InflationFactor,
                                 sfu=sfPower,sfupar=planned_bnds$gammaA,
                                 sfl=sfPower,sflpar=planned_bnds$gammaB)
    }
    
    details <- c(critical=StandardDesign$upper$bound[k])
    if(Z>StandardDesign$upper$bound[k]){
      decision="Efficacy"
      reason <- "Efficacy"
    } else {
      decision="Futility"
      reason <- "Futility"
    }
    
  }
  
  if(plot){
    if(analysis=="final"){
      warning("Cannot produce the plot yet, when analysis=final.")
    }else{
      PlotBoundaries(Bounds)
      points(analysis_res$Info/planned_bnds$Info.max,Z,pch=4,cex=1.5)  #at some point we may wish to add that all previous analyses are plotted as well
    }
  }
  
  out <- list(decision=decision,reason=reason,details=details)
  out
}

#----------------------------------------------------------------------
### Decision.R ends here
