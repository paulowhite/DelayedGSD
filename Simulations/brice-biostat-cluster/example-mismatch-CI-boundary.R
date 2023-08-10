### example-debug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 13 2023 (09:28) 
## Version: 
## Last-Updated: aug  9 2023 (11:03) 
##           By: Brice Ozenne
##     Update #: 23
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

method <- 3
binding <- FALSE
seed <- 579018813 ## 996631745
fixC <- TRUE
ar <- 10
hypo <- "power"

thets <- c(140, 140, 140)
nGSD <- c(549, 548, 549)

debug.plannedB <- CalcBoundaries(kMax = 2,  
                                 alpha = 0.025, 
                                 beta = 0.2,  
                                 InfoR.i = c(0.58, 1.00),  
                                 InfoR.d = c(0.68,1),  
                                 rho_alpha = 2,  
                                 rho_beta = 2,  
                                 method = method,  
                                 cNotBelowFixedc = fixC,
                                 bindingFutility= binding,
                                 delta = 0.6)

df.sim <- GenData(n = max(nGSD), 
                  N.fw = 2,
                  rand.block = c(1,1,0,0),
                  allsd = c(2.5,2.1,2.4),
                  mean0 = c(10,0,0),
                  delta = c(0,0.3,0.6),
                  ar = (0.86*2)*2*10,
                  cor.01.1 = -0.15,
                  cor.ij.1 = 0.68,
                  cor.0j.1 = -0.27,
                  seed = seed,
                  MissProb = matrix(c(0.04807692, 0.05769231, 0.00961538, 0.88461538), 
                                    nrow = 2, 
                                    ncol = 2, 
                                    dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing")) 
                                    ),
                  DigitsOutcome = 2,
                  TimeFactor = 14,
                  DigitsTime = 0)$d

## ** interim
SelectData(df.sim,t=thets[method])

debug.lmmI <- analyzeData(SelectData(df.sim,t=thets[method]),
                          ddf = "nlme", data.decision = sum(df.sim$t1 <= thets[method] + 1.50001*14), getinfo = TRUE, trace = TRUE)
debug.GSDI <- update(debug.plannedB, delta = debug.lmmI, trace = TRUE)
## summary(debug.GSDI)

##  Boundaries and observed statistics 
##  stage |         Interim           | Decision       |     Spent          
##        | F-bound E-bound    Stat   |  C-bound Stat  |     alpha      beta
##      1 |  0.9117 2.31555 2.17459 C |  1.95996       | 0.0102914 0.0823316
##      2 |                           |                |                    


## ** decision
debug.lmmD <- analyzeData(df.sim[which(df.sim$t1 <= thets[method] + 1.50001*14),],
                          ddf = "nlme", getinfo = TRUE, trace = TRUE)
debug.GSDD <- update(debug.GSDI, delta = debug.lmmD, k = 1, type.k = "decision", trace = FALSE)
summary(debug.GSDD)

##            GSD with repeated measurements at the decision analysis of stage 1 

## Boundaries and observed statistics 
## stage |         Interim             | Decision           | Spent     
##       | F-bound E-bound Stat        |  C-bound    Stat   | alpha beta
##     1 |    <NA>    <NA> <NA> S-Imax |  1.95996 1.73225 F | 0.025  0.2
##     2 |                             |                    |           

## Observed and predicted information: 
## stage |  Interim     (%) | Decision     (%) |   n
##     1 | 18.77808 0.81901 | 21.89924 0.95514 | 399
##     2 |                  |                  |    

## Current MUE-estimate of the treatment effect (Z1) 
## estimate   lower   upper  p.value
##   0.7975 0.47701 1.22493 1.000000

confint(debug.GSDD)

##----------------------------------------------------------------------
### example-debug.R ends here
