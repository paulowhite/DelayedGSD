### example-debug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 13 2023 (09:28) 
## Version: 
## Last-Updated: okt  6 2023 (17:37) 
##           By: Brice Ozenne
##     Update #: 70
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(DelayedGSD)
library(ggplot2)


## * Setting (for scenario 9)
binding <- FALSE
fixC <- TRUE
ar <- 10
hypo <- "power"
PropForInterim <- 0.5

nGSD <- c(557, 557, 557)

planned.M2 <- CalcBoundaries(kMax = 2,  
                             alpha = 0.025, 
                             beta = 0.2,  
                             InfoR.i = c(0.58, 1.00),  
                             InfoR.d = c(0.68,1),  
                             rho_alpha = 2,  
                             rho_beta = 2,  
                             method = 2,  
                             cNotBelowFixedc = fixC,
                             bindingFutility= binding,
                             delta = 0.6)
planned.M3 <- CalcBoundaries(kMax = 2,  
                             alpha = 0.025, 
                             beta = 0.2,  
                             InfoR.i = c(0.58, 1.00),  
                             InfoR.d = c(0.68,1),  
                             rho_alpha = 2,  
                             rho_beta = 2,  
                             method = 3,  
                             cNotBelowFixedc = fixC,
                             bindingFutility= binding,
                             delta = 0.6)

## * Example 1 (stopping for efficacy)

## ** generate data
seed1 <- 413883402 ##  413883402 or 341365344
df1.sim <- GenData(n = max(nGSD), 
                   N.fw = 2,
                   rand.block = c(1,1,0,0),
                   allsd = c(2.5,2.1,2.4),
                   mean0 = c(10,0,0),
                   delta = c(0,0.3,0.6),
                   ar = (0.86*2)*2*10,
                   cor.01.1 = -0.15,
                   cor.ij.1 = 0.68,
                   cor.0j.1 = -0.27,
                   seed = seed1,
                   MissProb = matrix(c(0.04807692, 0.05769231, 0.00961538, 0.88461538), 
                                     nrow = 2, 
                                     ncol = 2, 
                                     dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing")) 
                                     ),
                   DigitsOutcome = 2,
                   TimeFactor = 14,
                   DigitsTime = 0)$d
thets1 <- df1.sim$t3[ceiling(unique(nGSD)*PropForInterim)] ## df1.sim$t3[ceiling(nGSD*PropForInterim)]

## ** interim
debug1.lmmI <- analyzeData(SelectData(df1.sim,t=thets1),
                           ddf = "nlme", data.decision = sum(df1.sim$t1 <= thets1 + 1.50001*14), getinfo = TRUE, trace = TRUE)

debug1.IM2 <- update(planned.M2, delta = debug1.lmmI, trace = TRUE)
summary(debug1.IM2)
debug1.IM3 <- update(planned.M3, delta = debug1.lmmI, trace = TRUE)
summary(debug1.IM3)

##            GSD with repeated measurements at the interim analysis of stage 1 

## Boundaries and observed statistics 
## stage |         Interim             | Decision       |     Spent          
##       | F-bound E-bound    Stat     |  C-bound Stat  |     alpha      beta
##     1 | 0.48448 2.33168 2.70865 S-E |  1.95996       | 0.0073032 0.0584254
##     2 |                             |                |                    

## Observed and predicted information: 
## stage |  Interim     (%) | Decision    (%) |        n
##     1 | 12.58547 0.54049 | 16.58381 0.7122 | 407.0000
##     2 |                  |                 | 753.0235

## Current ML-estimate of the treatment effect (Z1) 
## estimate      se   lower   upper statistic  df   p.value
##  0.76352 0.28188 0.20981 1.31723   2.70865 543 0.0034844
                                                                   
## 407 patients in the study: 300 with at least one outcome value 
##                               249 with complete data              


## ** decision
debug1.lmmD <- analyzeData(df1.sim[which(df1.sim$t1 <= thets1 + 1.50001*14),],
                           ddf = "nlme", getinfo = TRUE, trace = TRUE)
debug1.DM2 <- update(debug1.IM2, delta = debug1.lmmD, k = 1, type.k = "decision", trace = FALSE,
                     continuity.correction = 2)
summary(debug1.DM2)
##  	    GSD with repeated measurements at the decision analysis of stage 1 

##  * Boundaries and observed statistics 
##  stage |         Interim             | Decision           |     Spent          
##        | F-bound E-bound    Stat     |  C-bound    Stat   |     alpha      beta
##      1 | 0.51244 2.44745 2.70865 S-E |  1.95996 1.89075 F | 0.0071936 0.0575486
##      2 |                             |                    |                    

##  * Observed and predicted information: 
##  stage |  Interim     (%) | Decision     (%) |   n
##      1 | 12.58547 0.53642 | 18.28058 0.77915 | 407
##      2 |                  |                  |    

##  * Current MUE-estimate of the treatment effect (Z1) 
##  estimate    lower   upper  p.value
##   0.90993 -0.11562 1.37777 1.000000
                                                                   
##   * 407 patients in the study: 389 with at least one outcome value 
##                                363 with complete data              

## Advarselsbesked:
## I confint.delayedGSD(x, method = c("ML", "MUE")) :
##   Possibly incorrect evaluation of the MUE lower, estimate.
## Mismatch in the optimization process in term of confidence level: 0.950416370875677, 0.25. 

debug1.DM3 <- update(debug1.IM3, delta = debug1.lmmD, k = 1, type.k = "decision", trace = FALSE,
                     continuity.correction = 2)
summary(debug1.DM3)

## 	    GSD with repeated measurements at the decision analysis of stage 1 

##  Boundaries and observed statistics 
##  stage |         Interim             | Decision           |     Spent          
##        | F-bound E-bound    Stat     |  C-bound    Stat   |     alpha      beta
##      1 | 0.48448 2.33168 2.70865 S-E |  1.95996 1.89075 F | 0.0073032 0.0584254
##      2 |                             |                    |                    

##  Observed and predicted information: 
##  stage |  Interim     (%) | Decision     (%) |   n
##      1 | 12.58547 0.54049 | 18.28058 0.78507 | 407
##      2 |                  |                  |    

##  Current MUE-estimate of the treatment effect (Z1) 
##  estimate  lower   upper  p.value
##    0.3047 0.3047 0.32528 0.997904
                                                                   
##   407 patients in the study: 389 with at least one outcome value 
##                              363 with complete data              

## Advarselsbesked:
## I confint.delayedGSD(x, method = c("ML", "MUE")) :
##   Possibly incorrect evaluation of the MUE lower, estimate.
## Mismatch in the optimization process in term of confidence level: 0.940666047062586, 0.244905618768506. 

## ** debug
calcP_new(delta = 0, object = debug1.DM2)
## [1] 0.00965248
## attr(,"error")
## [1] 1e-15
## attr(,"msg")
## [1] "Normal Completion"
## attr(,"terms")
## [1] 0.000000000 0.001466592 0.008185888

calcP_new(delta = 0, object = debug1.DM3)
## [1] 0.9979044
## attr(,"error")
## [1] 0
## attr(,"msg")
## [1] "univariate: using pnorm"
## attr(,"terms")
## [1] 0.9901412723 0.0004599626 0.0073031769

seqDelta <- seq(-5,5, length.out = 100)
dfW1.seq <- as.data.frame(do.call(rbind,lapply(seqDelta, function(iD){ ## iD <- 2
    iOut2 <- calcP_new(delta = iD, object = debug1.DM2)
    iOut3 <- calcP_new(delta = iD, object = debug1.DM3)

    iVec2 <- setNames(c(iD,2,iOut2,attr(iOut2,"terms")),
                     c("delta","method","p.value",paste0("term",1:3)))
    iVec3 <- setNames(c(iD,3,iOut3,attr(iOut3,"terms")),
                     c("delta","method","p.value",paste0("term",1:3)))

    return(rbind(iVec2,iVec3))
})))
dfL1.seq <- reshape(dfW1.seq, direction = "long", idvar = c("method","delta"), varying = list(c("p.value",paste0("term",1:3))),
                   timevar = "type", v.names = "value", times = c("p.value",paste0("term",1:3)))
rownames(dfL1.seq) <- NULL

ggP.ex1 <- ggplot(dfL1.seq, aes(x = delta, y = value, group = type, color = type, linetype = type)) 
ggP.ex1 <- ggP.ex1 + geom_line(linewidth=1.2) + geom_point(size=2) + facet_wrap(~method, labeller = label_both)
ggP.ex1 <- ggP.ex1 + scale_linetype_manual(values=c(1,2,2,2))
ggP.ex1

## * Example 2 (stopping for futility)

## ** generate data
seed2 <- 996631745 ##  996631745 or 657706742
df2.sim <- GenData(n = max(nGSD), 
                   N.fw = 2,
                   rand.block = c(1,1,0,0),
                   allsd = c(2.5,2.1,2.4),
                   mean0 = c(10,0,0),
                   delta = c(0,0.3,0.6),
                   ar = (0.86*2)*2*10,
                   cor.01.1 = -0.15,
                   cor.ij.1 = 0.68,
                   cor.0j.1 = -0.27,
                   seed = seed2,
                   MissProb = matrix(c(0.04807692, 0.05769231, 0.00961538, 0.88461538), 
                                     nrow = 2, 
                                     ncol = 2, 
                                     dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing")) 
                                     ),
                   DigitsOutcome = 2,
                   TimeFactor = 14,
                   DigitsTime = 0)$d
thets2 <- df2.sim$t3[ceiling(unique(nGSD)*PropForInterim)] ## df2.sim$t3[ceiling(nGSD*PropForInterim)]

## ** interim
debug2.lmmI <- analyzeData(SelectData(df2.sim,t=thets2),
                           ddf = "nlme", data.decision = sum(df2.sim$t1 <= thets2 + 1.50001*14), getinfo = TRUE, trace = TRUE)
debug2.IM2 <- update(planned.M2, delta = debug2.lmmI, trace = TRUE)
summary(debug2.IM2)
debug2.IM3 <- update(planned.M3, delta = debug2.lmmI, trace = TRUE)
summary(debug2.IM3)

##     GSD with repeated measurements at the interim analysis of stage 1 

##  * Boundaries and observed statistics 
##  stage |         Interim             | Decision       |     Spent          
##        | F-bound E-bound    Stat     |  C-bound Stat  |     alpha      beta
##      1 |   0.709 2.24163 0.59188 S-F |  1.95996       | 0.0090704 0.0725629
##      2 |                             |                |                    

##  * Observed and predicted information: 
##  stage |  Interim     (%) | Decision     (%) |        n
##      1 | 14.02574 0.60234 |  17.8922 0.76839 | 404.0000
##      2 |                  |                  | 670.7168

##  * Current ML-estimate of the treatment effect (Z1) 
##  estimate      se    lower   upper statistic  df p.value
##   0.15804 0.26702 -0.36647 0.68255   0.59188 543 0.27709
                                                                   
##   * 404 patients in the study: 308 with at least one outcome value 
##                                241 with complete data              


## ** decision
debug2.lmmD <- analyzeData(df2.sim[which(df2.sim$t1 <= thets2 + 1.50001*14),],
                           ddf = "nlme", getinfo = TRUE, trace = TRUE)
debug2.DM2 <- update(debug2.IM2, delta = debug2.lmmD, k = 1, type.k = "decision", trace = FALSE,
                     continuity.correction = 2)
summary(debug2.DM2)
## GSD with repeated measurements at the decision analysis of stage 1 

##  * Boundaries and observed statistics 
##  stage |         Interim             | Decision           |     Spent          
##        | F-bound E-bound    Stat     |  C-bound    Stat   |     alpha      beta
##      1 | 0.74485 2.36833 0.59188 S-F |  1.95996 2.33599 E | 0.0089342 0.0714739
##      2 |                             |                    |                    

##  * Observed and predicted information: 
##  stage |  Interim    (%) | Decision     (%) |   n
##      1 | 14.02574 0.5978 | 17.97582 0.76616 | 404
##      2 |                 |                  |    

##  * Current MUE-estimate of the treatment effect (Z1) 
##  estimate   lower   upper   p.value
##   0.63796 0.13268 1.15625 0.0064764
                                                                   
##   * 404 patients in the study: 382 with at least one outcome value 
##                                354 with complete data              

debug2.DM3 <- update(debug2.IM3, delta = debug2.lmmD, k = 1, type.k = "decision", trace = FALSE,
                     continuity.correction = 2)
summary(debug2.DM3)

## 	    GSD with repeated measurements at the decision analysis of stage 1 

##  Boundaries and observed statistics 
##  stage |         Interim             | Decision           |     Spent          
##        | F-bound E-bound    Stat     |  C-bound    Stat   |     alpha      beta
##      1 |   0.709 2.24163 0.59188 S-F |  1.95996 2.33599 F | 0.0090704 0.0725629
##      2 |                             |                    |                    

##  Observed and predicted information: 
##  stage |  Interim     (%) | Decision     (%) |   n
##      1 | 14.02574 0.60234 | 17.97582 0.77198 | 404
##      2 |                  |                  |    

##  Current MUE-estimate of the treatment effect (Z1) 
##  estimate   lower   upper   p.value
##   0.40601 0.40603 0.43304 0.9965833
                                                                   
##  404 patients in the study: 382 with at least one outcome value 
##                                354 with complete data              

## Advarselsbesked:
## I confint.delayedGSD(x, method = c("ML", "MUE")) :
##   Possibly incorrect evaluation of the MUE lower, estimate.
## Mismatch in the optimization process in term of confidence level: 0.923360387330875, 0.236114777922808. 

## ** debug
calcP_new(delta = 0, object = debug2.DM2)

calcP_new(delta = 0, object = debug2.DM3)
## [1] 0.9931424
## attr(,"error")
## [1] 0
## attr(,"msg")
## [1] "univariate: using pnorm"
## attr(,"terms")
## [1] 9.875072e-01 5.629446e-03 5.718868e-06


## [1] 0.9979044
## attr(,"error")
## [1] 0
## attr(,"msg")
## [1] "univariate: using pnorm"
## attr(,"terms")
## [1] 0.9901412723 0.0004599626 0.0073031769

seqDelta <- seq(-5,5, length.out = 100)
dfW2.seq <- as.data.frame(do.call(rbind,lapply(seqDelta, function(iD){ ## iD <- 2
    iOut2 <- calcP_futility(delta = iD, object = debug1.DM2)
    iOut3 <- calcP_futility(delta = iD, object = debug1.DM3)

    iVec2 <- setNames(c(iD,2,iOut2,attr(iOut2,"terms")),
                     c("delta","method","p.value",paste0("term",1:3)))
    iVec3 <- setNames(c(iD,3,iOut3,attr(iOut3,"terms")),
                     c("delta","method","p.value",paste0("term",1:3)))

    return(rbind(iVec2,iVec3))
})))
dfL2.seq <- reshape(dfW2.seq, direction = "long", idvar = c("method","delta"), varying = list(c("p.value",paste0("term",1:3))),
                   timevar = "type", v.names = "value", times = c("p.value",paste0("term",1:3)))
rownames(dfL2.seq) <- NULL

ggP.ex2 <- ggplot(dfL.seq, aes(x = delta, y = value, group = type, color = type, linetype = type)) 
ggP.ex2 <- ggP.ex2 + geom_line(linewidth=1.2) + geom_point(size=2) + facet_wrap(~method, labeller = label_both)
ggP.ex2 <- ggP.ex2 + scale_linetype_manual(values=c(1,2,2,2))
ggP.ex2
##----------------------------------------------------------------------
### example-debug.R ends here
