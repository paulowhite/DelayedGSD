### figure-pvalue-correction.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 10 2023 (13:16) 
## Version: 
## Last-Updated: mar 23 2023 (10:54) 
##           By: Brice Ozenne
##     Update #: 30
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(DelayedGSD)
library(RColorBrewer)
export <- TRUE

## * 2 stages

## ** graphical display - no fix C
if(export){
  pdf("figures/illustration-pvalue-2stage.pdf", width = 5)
}
res2stage <- gridFinalPvalue(list(Info.d = c(12.58797814, 20.74866926),
                                  Info.i = c(10.60271846),
                                  ck = c(1.49772992,1.9961649), ck.unrestricted = c(1.49772992,1.9961649),
                                  lk = c(0.24335814),#, 1.9961649),
                                  uk = c(2.5458844),#, 1.9961649),
                                  kMax = 2,
                                  reason.interim = c("no boundary crossed",NA),
                                  method = 3,
                                  bindingFutility = TRUE, 
                                  cNotBelowFixedc = FALSE),
                             continuity.correction = 0,
                             xlim = c(0.7,2.3))
if(export){
  dev.off()
}

## ** graphical display - fix C
res2stage.fixC <- vector(mode = "list", length = 3)

if(export){
  pdf("figures/illustration-pvalue-2stage-fixC.pdf", width = 12)
}
par(mfrow = c(1,3), mar = rep(4,4))
for(iC in 0:2){ ## iC <- 0
  iTitle <- switch(as.character(iC),
                   "0" = "no correction",
                   "1" = "continuity correction \n (add P)",
                   "2" = "continuity correction \n (shift stat)",
  )
  set.seed(10)
  res2stage.fixC[[iC+1]] <- gridFinalPvalue(list(Info.d = c(12.58797814),
                                                 Info.i = c(10.60271846, 20.74866926),
                                                 ck = c(1.959964,1.9961649), ck.unrestricted = c(1.49772992,1.9961649),
                                                 lk = c(0.24335814, 1.9961649),
                                                 uk = c(2.5458844, 1.9961649),
                                                 kMax = 2,
                                                 reason.interim = c("no boundary crossed",NA),
                                                 method = 3,
                                                 bindingFutility = TRUE, 
                                                 cNotBelowFixedc = TRUE),
                                            continuity.correction = iC,
                                            xlim = c(0.7,2.2), title = iTitle)
}
if(export){
  dev.off()
}

## * 3 stages

## ** graphical display - no fix C
bound3s <- CalcBoundaries(kMax = 3,  
                          alpha = 0.025, 
                          beta = 0.2,  
                          InfoR.i = c(0.40, 0.7, 1.00),  
                          InfoR.d = c(0.45, 0.75, 1),  
                          rho_alpha = 2,  
                          rho_beta = 2,  
                          method = 1,  
                          cNotBelowFixedc = FALSE,
                          bindingFutility= TRUE,
                          delta = 0.6)

## Planned boundaries: 
## stage  F-bound E-bound C-bound alpha-spent beta-spent
##     1 -0.01468 2.65207 1.39026     0.00400      0.032
##     2  1.07723 2.32431 1.74696     0.01225      0.098
##     3                  2.03284     0.02500      0.200

## Planned information: 
## stage  Interim (%) Decision  (%)
##     1  9.37888 0.4 10.55124 0.45
##     2 16.41305 0.7 17.58541 0.75
##     3              23.44721 1.00

if(export){
  pdf("figures/illustration-pvalue-3stage.pdf", width = 6)
}
res3stage <- gridFinalPvalue(bound3s,
                             continuity.correction = 0,
                             xlim = c(0.6,3.4))
if(export){
  dev.off()
}


## ** graphical display - fix C
bound3s.fixC <- CalcBoundaries(kMax = 3,  
                               alpha = 0.025, 
                               beta = 0.2,  
                               InfoR.i = c(0.40, 0.7, 1.00),  
                               InfoR.d = c(0.45, 0.75, 1),  
                               rho_alpha = 2,  
                               rho_beta = 2,  
                               method = 1,  
                               cNotBelowFixedc = TRUE,
                               bindingFutility= TRUE,
                               delta = 0.6)

res3stage.fixC <- vector(mode = "list", length = 3)
## Planned boundaries: 
## stage  F-bound E-bound C-bound alpha-spent beta-spent
##     1 -0.01468 2.65207 1.95996     0.00400      0.032
##     2  1.07723 2.32431 1.95996     0.01225      0.098
##     3                  2.03284     0.02500      0.200

## Planned information: 
## stage  Interim (%) Decision  (%)
##     1  9.37888 0.4 10.55124 0.45
##     2 16.41305 0.7 17.58541 0.75
##     3              23.44721 1.00

if(export){
  pdf("figures/illustration-pvalue-3stage-fixC.pdf", width = 12)
}
par(mfrow = c(1,3), mar = rep(4,4))
for(iC in 0:2){ ## iC <- 0
  iTitle <- switch(as.character(iC),
                   "0" = "no correction",
                   "1" = "continuity correction \n (add P)",
                   "2" = "continuity correction \n (shift stat)",
  )
  set.seed(10)
  res3stage.fixC[[iC+1]] <- gridFinalPvalue(bound3s.fixC,
                                            continuity.correction = iC,
                                            xlim = c(0.65,3.3),
                                            title = iTitle,
                                            digits = 3)
}
if(export){
  dev.off()
}



##----------------------------------------------------------------------
### figure-pvalue-correction.R ends here