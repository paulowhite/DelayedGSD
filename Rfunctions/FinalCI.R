## * FinalCI
FinalCI <- function(Id,  #Information at all decision analyses up to stage where study was stopped
                    Ik,  #Information at all interim analyses up to stage where study was stopped
                    ck,   #decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis)
                    lk,  #lower bounds up to stage where study was stopped
                    uk,  #upper bounds up to stage where study was stopped
                    sided=1,  #whether the test is 1 or 2-sided
                    kMax, #maximum number of analyses
                    alpha=0.05, #confidence level (to get a 100*(1-alpha)% CI)
                    estimate,
                    optimizer = "optimise"){ #the observed treatment estimate at decision

  f <- function(delta){
    FinalPvalue(Id=Id,
                 Ik=Ik,
                 ck=ck,
                 lk=lk,
                 uk=uk,
                 sided=sided,
                 kMax=kMax,
                 estimate=estimate,
                 delta=delta) - alpha/2
  }

  g <- function(delta){
    1-FinalPvalue(Id=Id,
                 Ik=Ik,
                 ck=ck,
                 lk=lk,
                 uk=uk,
                 sided=sided,
                 kMax=kMax,
                 estimate=estimate,
                 delta=delta) - alpha/2
  }

    if(optimizer=="optimise"){
        lbnd <- optimise(function(x){f(x)^2},lower=estimate-4*sqrt(1/Id[length(Id)]),upper=estimate+0.5*sqrt(1/Id[length(Id)]))
        ubnd <- optimise(function(x){g(x)^2},lower=estimate-0.5*sqrt(1/Id[length(Id)]),upper=estimate+4*sqrt(1/Id[length(Id)]))
        out <- c(lower = lbnd$minimum, upper = ubnd$minimum)
        attr(out,"error") <- c(lower = lbnd$objective, upper = ubnd$objective)
    }else if(optimizer=="uniroot"){
        lbnd <- uniroot(f,lower=estimate-4*sqrt(1/Id[length(Id)]),upper=estimate+0.5*sqrt(1/Id[length(Id)]))
        ubnd <- uniroot(g,lower=estimate-0.5*sqrt(1/Id[length(Id)]),upper=estimate+4*sqrt(1/Id[length(Id)]))
        out <- c(lower = lbnd$root, upper = ubnd$root)
        attr(out,"error") <- c(lower = lbnd$f.root, upper = ubnd$f.root)
        attr(out,"iter") <- c(lower = lbnd$iter, upper = ubnd$iter)
    }

  return(out)
}
