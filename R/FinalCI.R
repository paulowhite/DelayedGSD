## * FinalCI (documentation)
#' @title calculate confidence intervals at the end of the study
#' 
#' @param Info.d Information at all decision analyses up to stage where study was stopped
#' @param Info.i Information at all interim analyses up to stage where study was stopped
#' @param ck decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis)
#' @param lk lower bounds up to stage where study was stopped
#' @param uk upper bounds up to stage where study was stopped
#' @param kMax maximum number of analyses
#' @param alpha confidence level (to get a 100*(1-alpha)\% CI)
#' @param estimate naive estimate (e.g. using  ML or REML).
#' @param optimizer the observed treatment estimate at decision
#' @param method  method 1, 2 or 3
#' @param bindingFutility [logical]  whether the futility stopping rule is binding.
#' @param cNotBelowFixedc [logical] whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)

## * FinalCI (code)
#' @export
FinalCI <- function(Info.d,  
                    Info.i,  
                    ck,  
                    lk,  
                    uk,  
                    kMax, 
                    alpha=0.05,
                    estimate,
                    optimizer = "optimise",
                    method,
                    bindingFutility,
                    cNotBelowFixedc){

    f <- function(delta){
        FinalPvalue(Info.d=Info.d,
                    Info.i=Info.i,
                    ck=ck,
                    lk=lk,
                    uk=uk,
                    kMax=kMax,
                    estimate=estimate,
                    delta=delta,
                    method=method,
                    bindingFutility=bindingFutility,
                    cNotBelowFixedc=cNotBelowFixedc) - alpha/2
    }

    g <- function(delta){
        1-FinalPvalue(Info.d=Info.d,
                      Info.i=Info.i,
                      ck=ck,
                      lk=lk,
                      uk=uk,
                      kMax=kMax,
                      estimate=estimate,
                      delta=delta,
                      method=method,
                      bindingFutility=bindingFutility,
                      cNotBelowFixedc=cNotBelowFixedc) - alpha/2
    }

    if(optimizer=="optimise"){
        lbnd <- stats::optimise(function(x){f(x)^2},lower=estimate-4*sqrt(1/Info.d[length(Info.d)]),upper=estimate+0.5*sqrt(1/Info.d[length(Info.d)]), tol = 1e-10)
        ubnd <- stats::optimise(function(x){g(x)^2},lower=estimate-0.5*sqrt(1/Info.d[length(Info.d)]),upper=estimate+4*sqrt(1/Info.d[length(Info.d)]), tol = 1e-10)
        out <- c(lower = lbnd$minimum, upper = ubnd$minimum)
        attr(out,"error") <- c(lower = lbnd$objective, upper = ubnd$objective)
    }else if(optimizer=="uniroot"){
        lbnd <- stats::uniroot(f,lower=estimate-4*sqrt(1/Info.d[length(Info.d)]),upper=estimate+0.5*sqrt(1/Info.d[length(Info.d)]), tol = 1e-10)
        ubnd <- stats::uniroot(g,lower=estimate-0.5*sqrt(1/Info.d[length(Info.d)]),upper=estimate+4*sqrt(1/Info.d[length(Info.d)]), tol = 1e-10)
        out <- c(lower = lbnd$root, upper = ubnd$root)
        attr(out,"error") <- c(lower = lbnd$f.root, upper = ubnd$f.root)
        attr(out,"iter") <- c(lower = lbnd$iter, upper = ubnd$iter)
    }

  return(out)
}
