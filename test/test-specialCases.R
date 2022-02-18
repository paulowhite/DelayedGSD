### test-specialCases.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  7 2022 (15:00) 
## Version: 
## Last-Updated: feb 18 2022 (12:24) 
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

library(testthat)

## * Skipped interim analysis
GSD.3t <- CalcBoundaries(kMax=3,  ## max number of analyses (including final)
                         alpha=0.025,  ## type I error
                         beta=0.2,  ## type II error
                         InfoR.i=c(0.5,0.75,1),  ## planned or observed information rates
                         rho_alpha=2,  ## rho parameter for alpha error spending function
                         rho_beta=2,  ## rho parameter for beta error spending function
                         method=1,  ## use method 1 or 2 from paper H&J
                         delta=1.5,  ## effect that the study is powered for
                         InfoR.d=c(0.55,0.8))

set.seed(1322)
GD <- GenData(n=82*2,delta=1.3,ar=5)  #generate data with all data for in case trial completes

lmmI1.3t <- analyzeData(GD$d[1:50,], ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDI1.3t <- update(GSD.3t, delta = lmmI1.3t)


update(GSD.3t, delta = lmmI1.3t)
updateBoundaries(GSD.3t, delta = 1.6777011, Info.i = 1.11657, Info.d = 1.11657, k = 1, type.k = "interim", update.stage = TRUE)
updateMethod1(rho_alpha = 2,
              rho_beta = 2,
              alpha = 0.025, alphaSpent = rep(NA,3),
              beta = 0.2, betaSpent = rep(NA,3),
              Kmax = 3,
              Info.max = 3.773015,
              uk = rep(NA,3),
              lk = rep(NA,3),
              k = 1, type.k = "interim", ImaxAnticipated = FALSE,
              InfoR.i = c(0.2959368, NA, NA),
              InfoR.d = c(0.2959368, NA, NA),
              delta = 1.6777011, 
              alternative = "greater",
              binding = TRUE,
              Trace = 0,
              cMin = -Inf)


lmmI2.3t <- analyzeData(GD$d[1:45,], ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDI2.3t <- update(GSDI1.3t, delta = lmmI2.3t)

plot(GSD.3t)
plot(GSDI1.3t)
plot(GSDI2.3t)
plot(GSDI2.3t, cex.estimate = c(2,0.7))



##----------------------------------------------------------------------
### test-specialCases.R ends here
