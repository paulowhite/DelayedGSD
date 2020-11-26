### Decision.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 28 2020 (14:40) 
## Version: 
## Last-Updated: Aug 28 2020 (15:28) 
##           By: Paul Blanche
##     Update #: 18
#----------------------------------------------------------------------
## 
### Commentary: 
##  Fucntion that maps the statistical results at interim (or final) analysis into
##  decision to reject, or continue or stop inclusind subjects.
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


##' @title xx
##' @description yy
##' @param object of class ??
##' @param rho paramter of the error spending functions
##' @param alpha type-I error of the two sided test
##' @details zz
##' @return ff
##' @author Paul Blanche
##' 
##' @examples
##' x <- rnorm(10)
##' @export
Decision <- function(analysis_res,  #results from AnalyzeData
                     planned_bnds,  #results from CalcBoundaries
                     k=1, #at which phase are we?
                     analysis="interim", #is it an interim or decision or final analysis
                     Ik,  #all I_k (information) from first interim to final analysis (observed where possible)
                     Id,  #expected or observed information at each decision analysis
                     plot=T){  #should the boundaries and results be plotted?
  
  Bounds <- CalcBoundaries(kMax=planned_bnds$kMax,
                           sided=planned_bnds$sided,  #one or two-sided
                           alpha=planned_bnds$alpha,  #type I error
                           beta=planned_bnds$beta,  #type II error
                           informationRates=Ik/planned_bnds$Imax,  #planned or observed information rates
                           gammaA=planned_bnds$gammaA,  #rho parameter for alpha error spending function
                           gammaB=planned_bnds$gammaB,  #rho parameter for beta error spending function
                           method=planned_bnds$method,  #use method 1 or 2 from paper H&J
                           delta=planned_bnds$delta,  #effect that the study is powered for
                           Id=Id)  #(expected) information ratio at each decision analysis
  
  Z <- analysis_res$estimate/analysis_res$se
  
  if(analysis=="interim"){
    if(Z>Bounds$uk[k]){
      decision="Efficacy"
    } else if (Z<Bounds$lk[k]){
      decision="Futility"
    } else {
      decision="Continue"
    }
  } else if(analysis=="decision"){
    if(Z>Bounds$ck[k]){
      decision="Efficacy"
    } else {
      decision="Futility"
    }
  } else if(analysis=="final"){
    if(Z>Bounds$uk[k]){
      decision="Efficacy"
    } else {
      decision="Futility"
    }
  } else {
    stop("Please specify interim or decision of final for analysis type")
  }
  
  if(plot){
    PlotBoundaries(Bounds)
    points(analysis_res$Info/planned_bnds$Imax,Z)  #at some point we may wish to add that all previous analyses are plotted as well
  }
  
  decision
  
}

#would be nice if we can have a predict_info function to avoid that Id needs to be an argument



#----------------------------------------------------------------------
### Decision.R ends here
