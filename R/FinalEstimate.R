## * FinalEstimate (documentation)
#' @title calculate de-biased estimate
#'
#' @param Info.d Information at all decision analyses up to stage where study was stopped (should include information at final analysis if stopped at final analysis)
#' @param Info.i Information at all interim analyses up to stage where study was stopped
#' @param ck decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis)
#' @param lk lower bounds up to stage where study was stopped
#' @param uk upper bounds up to stage where study was stopped
#' @param sided one or two sided
#' @param kMax maximum number of analyses
#' @param estimate

## * FinalEstimate (code)
#' @export
FinalEstimate <- function(Info.d,  
                          Info.i,  
                          ck,   
                          lk,  
                          uk,  
                          sided=1,  
                          kMax, 
                          estimate){ 
  
  f <- function(delta){
    (FinalPvalue(Info.d=Info.d,
                Info.i=Info.i,
                ck=ck,
                lk=lk,
                uk=uk,
                sided=sided,
                kMax=kMax,
                estimate=estimate,
                delta=delta) - 0.5)^2
  }
  
  res <- optimise(f,lower=estimate-2*sqrt(1/Info.d[length(Info.d)]),upper=estimate+2*sqrt(1/Info.d[length(Info.d)]))
  
  res$minimum
}
