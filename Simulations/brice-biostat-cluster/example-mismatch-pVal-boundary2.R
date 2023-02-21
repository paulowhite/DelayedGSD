### example-mismatch-pVal-boundary2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 20 2023 (17:36) 
## Version: 
## Last-Updated: jan 27 2023 (14:14) 
##           By: Brice Ozenne
##     Update #: 6
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * example where concluding for futility but with significant p-value
CalcBoundaries(kMax = 2,
               alpha = 0.025,
               beta = 0.2,
               InfoR.i = c(0.58, 1),
               InfoR.d = c(0.68, 1),
               rho_alpha = 2,
               rho_beta = 2,
               method = 1, 
               cNotBelowFixedc = TRUE,
               bindingFutility = TRUE,
               delta = 0.6)

d <- GenData(n = 553,
             N.fw = 2,
             rand.block = c(1, 1, 0, 0),
             allsd = c(2.5, 2.1, 2.4),
             mean0 = c(10, 0, 0),
             delta = c(0, 0.3, 0.6),
             ar = 34.4, 
             cor.01.1 = -0.15,
             cor.ij.1 = 0.68,
             cor.0j.1 = -0.27,
             seed = 490137693, 
             MissProb = structure(c(0.0480769230769231, 0.0576923076923077, 
                                    0.00961538461538462, 0.884615384615385),
                                  dim = c(2L, 2L), dimnames = list(c("V1 missing", "V1 not missing"), c("V2 missing", "V2 not missing"))),
             DigitsOutcome = 2,
             TimeFactor = 14,
             DigitsTime = 0)$d

myBound0 <- CalcBoundaries(kMax=theK,
                           alpha=theAlpha,
                           beta=theBeta,
                           InfoR.i=c(0.5,1),
                           rho_alpha=2,
                           rho_beta=2,
                           method=1,
                           delta=theDelta,
                           InfoR.d=0.55,
                           cNotBelowFixedc=TRUE)
myBound0
    
theAR <- 10  #accrual rate (pt per month)
theDelay <- 0.7500001  #time in months to process data
tau.i <- theData$t3[theN + ceiling(theAR*theDelay)] #time point at which to do IA

#### Analyse data at the first interim ####
theInterimData <- SelectData(theData, t = 0.5*tau.i)
myLmmI <- analyzeData(theInterimData)
myInterim1 <- update(myBound0, delta = myLmmI, trace = FALSE) ## k = 1, analysis = "interim"
myInterim1

#### Analyse data at the final stage ####
theFinalData <- theData[theData$id %in% c(theInterimData$id,max(theInterimData$id)+1:2),]
myLmmF <- analyzeData(theFinalData)
myFinal <- update(myInterim1, delta = myLmmF, trace = FALSE) ## k = 2, analysis = "final"

## decreasing information
myLmmF$information - myLmmI$information

## issue with the calculation in updateMethod1:

sigmaZk.debug <- matrix(c(1,sqrt(myInterim1$Info.i[1]/myLmmF$information[1]),sqrt(myInterim1$Info.i[1]/myLmmF$information[1]),1),2,2)
sigmaZk.debug

uniroot(function(x){pmvnorm(lower = c(1.336678,x),
                            upper = c(2.172381,Inf),
                            mean=rep(0,2),
                            sigma= sigmaZk.debug,
                            abseps = 1e-06) - 0.01008653},
        lower = 1.336678,
        upper = 2.172381,
        tol = 1e-06)$root


## * how to get there
## xxx <- res2stage[decision=="futility" & p.value_MUE<0.025,
##                  .(scenario,missing,binding,fixC,ar,hypo,method,stage,type,statistic,p.value_ML,p.value_MUE,
##                    ck,decision,reason,seed)]
## xxx[,.N,by=c("seed","scenario")][which.max(N)]
## res2stage[seed == 490137693 & scenario == 5, .SD, .SDcols = names(xxx)]
## ##    scenario missing binding fixC ar  hypo method stage     type statistic p.value_ML p.value_MUE       ck decision   reason      seed
## ## 1:        5    TRUE    TRUE TRUE 10 power      1     1  interim  2.542933         NA          NA       NA     stop efficacy 490137693
## ## 2:        5    TRUE    TRUE TRUE 10 power      1     1 decision  1.918589 0.02771629 0.005601226 1.959964 futility     <NA> 490137693
## ## 3:        5    TRUE    TRUE TRUE 10 power      1     2    final        NA         NA          NA       NA     <NA>     <NA> 490137693
## ## 4:        5    TRUE    TRUE TRUE 10 power      2     1  interim  2.542933         NA          NA       NA     stop efficacy 490137693
## ## 5:        5    TRUE    TRUE TRUE 10 power      2     1 decision  1.918589 0.02771629 0.005502683 1.959964 futility     <NA> 490137693
## ## 6:        5    TRUE    TRUE TRUE 10 power      2     2    final        NA         NA          NA       NA     <NA>     <NA> 490137693
## ## 7:        5    TRUE    TRUE TRUE 10 power      3     1  interim  2.542933         NA          NA       NA     stop efficacy 490137693
## ## 8:        5    TRUE    TRUE TRUE 10 power      3     1 decision  1.918589 0.02771629 0.006722260 1.959964 futility     <NA> 490137693
## ## 9:        5    TRUE    TRUE TRUE 10 power      3     2    final        NA         NA          NA       NA     <NA>     <NA> 490137693

##----------------------------------------------------------------------
### example-mismatch-pVal-boundary2.R ends here
