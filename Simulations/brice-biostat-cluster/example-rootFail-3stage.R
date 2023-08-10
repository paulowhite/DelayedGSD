### example-rootFail-3stage.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  9 2023 (11:06) 
## Version: 
## Last-Updated: aug  9 2023 (11:23) 
##           By: Brice Ozenne
##     Update #: 11
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Simulate data (special case)
library(DelayedGSD)

method <- 2
binding <- TRUE
seed <- 579018813 ## 996631745
fixC <- FALSE
ar <- 34.4
hypo <- "typeI"

thets <- cbind("time.interim1" = c(112, 112, 112), 
               "time.interim2" = c(162, 162, 162)
               )
nGSD <- c(552, 552, 553)

debug.plannedB <- CalcBoundaries(kMax = 3,  
                                 alpha = 0.025, 
                                 beta = 0.2,  
                                 InfoR.i = c(0.4,0.65, 1.00),  
                                 InfoR.d = c(0.5,0.75,1),  
                                 rho_alpha = 2,  
                                 rho_beta = 2,  
                                 method = method,  
                                 cNotBelowFixedc = fixC,
                                 bindingFutility= binding,
                                 delta = 0.6)
## many warnings


df.sim <- GenData(n = max(nGSD), 
                  N.fw = 3,
                  rand.block = c(1,1,0,0),
                  allsd = c(2.5, 2.1, 2.4, 2.7),
                  mean0 = c(10, 0, 0, 0),
                  delta = c(0, 0, 0, 0),
                  ar = ar,
                  cor.01.1 = -0.15,
                  cor.ij.1 = 0.68,
                  cor.0j.1 = -0.27,
                  seed = seed,
                  MissProb = structure(c(0.00961538461538462, 0.0192307692307692, 0.0288461538461538, 
                                         0.0384615384615385, 0.00961538461538462, 0.0192307692307692, 
                                         0.00961538461538462, 0.865384615384615), dim = c(2L, 2L, 2L),
                                       dimnames = list(c("V1 missing", "V1 not missing"),
                                                       c("V2 missing", "V2 not missing"),
                                                       c("V3 missing", "V3 not missing"))),
                  DigitsOutcome = 2,
                  TimeFactor = 14,
                  DigitsTime = 0)$d


##----------------------------------------------------------------------
### example-rootFail-3stage.R ends here
