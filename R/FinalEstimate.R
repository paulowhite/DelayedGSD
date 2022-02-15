## * FinalEstimate (documentation)
#' @title calculate de-biased estimate
#'
#' @param Info.d Information at all decision analyses up to stage where study was stopped (should include information at final analysis if stopped at final analysis)
#' @param Info.i Information at all interim analyses up to stage where study was stopped
#' @param ck decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis)
#' @param lk lower bounds up to stage where study was stopped
#' @param uk upper bounds up to stage where study was stopped
#' @param kMax maximum number of analyses
#' @param estimate naive estimate (e.g. using  ML or REML).
#' @param futility2efficacy [logical] is it possible to stop for futility at interim and reject the null hypothesis at decision.
#' @param bindingFutility [logical]  whether the futility stopping rule is binding.

## * FinalEstimate (code)
#' @export
FinalEstimate <- function(Info.d,  
                          Info.i,  
                          ck,   
                          lk,  
                          uk,  
                          kMax, 
                          estimate,
                          futility2efficacy,
                          bindingFutility){ 
  
  f <- function(delta){
    (FinalPvalue(Info.d=Info.d,
                Info.i=Info.i,
                ck=ck,
                lk=lk,
                uk=uk,
                kMax=kMax,
                estimate=estimate,
                delta=delta,
                futility2efficacy=futility2efficacy,
                bindingFutility=bindingFutility) - 0.5)^2
  }
  
  res <- stats::optimise(f,lower=estimate-2*sqrt(1/Info.d[length(Info.d)]),upper=estimate+2*sqrt(1/Info.d[length(Info.d)]))
  
  res$minimum
}
