### test-specialCases.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  7 2022 (15:00) 
## Version: 
## Last-Updated: feb 23 2022 (14:16) 
##           By: Brice Ozenne
##     Update #: 10
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
GD <- GenData(n=82*2,delta=1.3,ar=5)  ## generate data with all data for in case trial completes
GD.I <- SelectData(GD$d, t = GD$d[50,"t3"])
    
lmmI1.3t <- analyzeData(GD.I, ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDI1.3t <- update(GSD.3t, delta = lmmI1.3t)

GD.I2 <- SelectData(GD$d, t = GD$d[45,"t3"])
lmmI2.3t <- analyzeData(GD.I2, ddf = "nlme", getinfo = TRUE, trace = TRUE)
GSDI2.3t <- update(GSDI1.3t, delta = lmmI2.3t)

plot(GSDI2.3t)
plot(GSDI2.3t, cex.estimate = c(2,0.7))



##----------------------------------------------------------------------
### test-specialCases.R ends here
