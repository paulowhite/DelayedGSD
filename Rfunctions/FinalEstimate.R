FinalEstimate <- function(Info.d,  #Information at all decision analyses up to stage where study was stopped (should include information at final analysis if stopped at final analysis)
                          Info.i,  #Information at all interim analyses up to stage where study was stopped
                          ck,   #decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis)
                          lk,  #lower bounds up to stage where study was stopped
                          uk,  #upper bounds up to stage where study was stopped
                          sided=1,  #whether the test is 1 or 2-sided
                          kMax, #maximum number of analyses
                          estimate){ #the observed treatment estimate at decision
  
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
