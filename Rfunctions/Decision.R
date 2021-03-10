### Decision.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 28 2020 (14:40) 
## Version: 
## Last-Updated: Mar  9 2021 (13:55) 
##           By: Paul Blanche
##     Update #: 69
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
                     Id=NULL,  #expected or observed information RATIO at each decision analysis
                     PositiveIsGood=TRUE, # whether positive effect is good (i.e. positive trial)
                     plot=T){  #should the boundaries and results be plotted?
    if(is.null(Id) & !((analysis=="interim" & planned_bnds$method==1) | analysis=="final")){
        stop("Id must be specified.")
    }
    if(is.null(Id) & analysis=="interim" & planned_bnds$method==1 & k==1){
        Id <- (Ik[1]/planned_bnds$Imax + 1)/2 #Rk: in this case it does not matter, this value is not used.
        ## print(paste0("Id=",Id))
    }

    if(analysis!="final"){
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
    }else{
        StandardDesign <- getDesignGroupSequential(kMax=planned_bnds$kMax,
                                                   sided = planned_bnds$sided,
                                                   alpha = planned_bnds$alpha,
                                                   beta = planned_bnds$beta,
                                                   informationRates = Ik/max(Ik),
                                                   typeOfDesign="asKD",
                                                   typeBetaSpending="bsKD",
                                                   gammaA=gammaA,
                                                   gammaB=gammaB)    
    }

    Z <- analysis_res$estimate/analysis_res$se
    if(!PositiveIsGood){
        Z <- -Z
        message("Negative effect is good: schange sign of Z")
    }
    if(analysis=="interim"){
        if(Z>Bounds$uk[k]){
            decision="Efficacy"
        } else if (Z<Bounds$lk[k]){
            decision="Futility"
        } else {
            decision="Continue"
        }
    } else if(analysis=="decision"){
        ## print(paste0("c=",Bounds$ck[k]))
        if(Z>Bounds$ck[k]){
            decision="Efficacy"
        } else {
            decision="Futility"
        }
    } else if(analysis=="final"){
        if(Z>StandardDesign$criticalValues[k]){
            decision="Efficacy"
        } else {
            decision="Futility"
        }
    } else {
        stop("Please specify interim or decision of final for analysis type")
    }
  
    if(plot){
        if(analysis=="final"){
            stop("Cannot produce the plot yet, when analysis=final.")
        }else{
            PlotBoundaries(Bounds)
            points(analysis_res$Info/planned_bnds$Imax,Z)  #at some point we may wish to add that all previous analyses are plotted as well
        }
    }
  
  decision
  
}

#would be nice if we can have a predict_info function to avoid that Id needs to be an argument



#----------------------------------------------------------------------
### Decision.R ends here
