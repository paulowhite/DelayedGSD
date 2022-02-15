### exampleGSD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Dec 16 2021 (08:56) 
## Version: 
## Last-Updated: Jan 28 2022 (16:00) 
##           By: Brice Ozenne
##     Update #: 29
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(DelayedGSD)

## * Design
thear <- 10  ## accrual rate (pt per month)
theDelta.t <- 0.7500001  ## time in months to process data

## planned boundaries
theAlpha <- 0.025
theBeta <- 0.2
theDelta <- 1.5

## * Boundaries
## ** 2 interim
GSD.2t <- CalcBoundaries(kMax=2,  ## max number of analyses (including final)
                         alpha=theAlpha,  ## type I error
                         beta=theBeta,  ## type II error
                         InfoR.i=c(0.5,1),  ## planned or observed information rates
                         rho_alpha=2,  ## rho parameter for alpha error spending function
                         rho_beta=2,  ## rho parameter for beta error spending function
                         method=2,  ## use method 1 or 2 from paper H&J
                         delta=theDelta,  ## effect that the study is powered for
                         InfoR.d=0.55,
                         cNotBelowFixedc=TRUE)

GSD.2t
coef(GSD.2t, type = "information")
coef(GSD.2t, type = "boundary")
print(GSD.2t)
plot(GSD.2t)

## ** 4 interim
GSD.4t <- CalcBoundaries(kMax=4,  ## max number of analyses (including final)
                         alpha=theAlpha,  ## type I error
                         beta=theBeta,  ## type II error
                         InfoR.i=c(0.25,0.5,0.75,1),  ## planned or observed information rates
                         rho_alpha=2,  ## rho parameter for alpha error spending function
                         rho_beta=2,  ## rho parameter for beta error spending function
                         method=1,  ## use method 1 or 2 from paper H&J
                         delta=theDelta,  ## effect that the study is powered for
                         InfoR.d=c(0.3,0.55,0.8))

GSD.4t
coef(GSD.4t, type = "information")
coef(GSD.4t, type = "boundary")
print(GSD.4t)
plot(GSD.4t)

## * Data
## ** simulation
set.seed(1322)
GD <- GenData(n=82*2,delta=1.3,ar=5)  #generate data with all data for in case trial completes

## ** time point at which to do the interim analysis
thet.2t <- GD$d$t3[ceiling(nrow(GD$d)/2) + ceiling(thear*theDelta.t)] ## 2 interim
## mean(x$d$t3<=thet.2t)
dataI.2t <- SelectData(GD$d,t=thet.2t)
dataF.2t <- GD$d[which(GD$t1 <= thet + theDelta.t*TimeFactor),]

thet.4t <- c(GD$d$t3[ceiling(nrow(GD$d)/4) + ceiling(thear*theDelta.t)],
             GD$d$t3[ceiling(nrow(GD$d)/2) + ceiling(thear*theDelta.t)], 
             GD$d$t3[ceiling(3*nrow(GD$d)/4) + ceiling(thear*theDelta.t)]) ## 4 interim,''
## c(mean(GD$d$t3<=thet.4t[1]),mean(GD$d$t3<=thet.4t[2]),mean(GD$d$t3<=thet.4t[3]))
dataI1.4t <- SelectData(GD$d,t=thet.4t[1])
dataI2.4t <- SelectData(GD$d,t=thet.4t[2])
dataI3.4t <- SelectData(GD$d,t=thet.4t[3])
dataF.4t <- GD$d

## * First interim
## ** 2 interim
GSDI.2t <- update(GSD.2t, data = dataI.2t)
print(GSDI.2t, planned = FALSE)
## coef(GSDI.2t, type = "effect")
## coef(GSDI.2t, type = "information", planned = FALSE)
## coef(GSDI.2t, type = "information", planned = TRUE)
## coef(GSDI.2t, type = "information", planned = "only")
## coef(GSDI.2t, type = "boundary", planned = FALSE)
## coef(GSDI.2t, type = "boundary", planned = TRUE)
plot(GSDI.2t)

GSDF.2t <- update(GSDI.2t, data = dataF.2t)

## ** 4 interim (decreasing)
GSDI1.4t <- update(GSD.4t, data = dataI1.4t)
GSDI1.4t$conclusion
GSDI1.4t$alphaSpent
GSDI1.4t$planned$alphaSpent

GSDI2.4t <- update(GSDI1.4t, data = dataI1.4t[-(1:10),])
GSDI2.4t$alphaSpent
GSDI2.4t$planned$alphaSpent

GSDI2.4t$conclusion
GSDI2.4t$delta
plot(GSDI2.4t)

GSDI3.4t <- update(GSDI2.4t, data = dataI3.4t)

GSDI3D.4t <- update(GSDI3.4t, data = dataI3.4t[-(1:20),])



GSDI3.4t$alphaSpent
GSDI3.4t$planned$alphaSpent

## ** 4 interim
GSDI1.4t <- update(GSD.4t, data = dataI1.4t)
print(GSDI1.4t)
print(GSDI1.4t, planned = TRUE)
## coef(GSDI1.4t, type = "effect")
## coef(GSDI1.4t, type = "information", planned = FALSE)
## coef(GSDI1.4t, type = "information", planned = TRUE)
## coef(GSDI1.4t, type = "information", planned = "only")
## coef(GSDI1.4t, type = "boundary", planned = FALSE)
## coef(GSDI1.4t, type = "boundary", planned = TRUE)
plot(GSDI1.4t)

## * Second interim
GSDI2.4t <- update(GSDI1.4t, data = dataI2.4t)
print(GSDI2.4t)
## coef(GSDI2.4t, type = "effect")
## coef(GSDI2.4t, type = "information", planned = FALSE)
## coef(GSDI2.4t, type = "information", planned = TRUE)
## coef(GSDI2.4t, type = "information", planned = "only")
## coef(GSDI2.4t, type = "boundary", planned = FALSE)
## coef(GSDI2.4t, type = "boundary", planned = TRUE)
## coef(GSDI2.4t, type = "decision")
plot(GSDI2.4t)

## * Third interim
## ISSUE: predicted boundary at decision is NA because too large information!
GSDI3.4t <- update(GSDI2.4t, data = dataI3.4t)
## coef(GSDI3.4t, type = "effect")
## coef(GSDI3.4t, type = "information")
## coef(GSDI3.4t, type = "boundary")
print(GSDI3.4t) 
plot(GSDI3.4t, planned = "only")
plot(GSDI3.4t, planned = TRUE)
plot(GSDI3.4t, planned = FALSE)



##----------------------------------------------------------------------
### exampleGSD.R ends here
