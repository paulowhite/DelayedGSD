### example-debug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 13 2023 (09:28) 
## Version: 
## Last-Updated: jan 13 2023 (15:04) 
##           By: Brice Ozenne
##     Update #: 13
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

theAlpha <- 0.025
theBeta <- 0.2
theDelta <- 1.5
theK <- 2
theN <- 82

theData <- GenData(n=theN*2,delta=theDelta*0,ar=5, seed = 76)$d
 
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

##----------------------------------------------------------------------
### example-debug.R ends here
