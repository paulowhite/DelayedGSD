### example-debug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 13 2023 (09:28) 
## Version: 
## Last-Updated: aug 10 2023 (18:32) 
##           By: Brice Ozenne
##     Update #: 31
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##    scenario method     type stage statistic     info       uk        lk       ck decision
## 1:        9      3  interim     1 0.5918775 14.02574 2.241626 0.7090032       NA     stop
## 2:        9      3 decision     1 2.3359933 17.97582       NA        NA 1.959964 futility
## 3:        9      3  interim     1 2.4761300 11.72226 2.416131 0.3589753       NA     stop
## 4:        9      3 decision     1 1.9329332 16.58190       NA        NA 1.959964 futility
##                          reason p.value_MUE lower_MUE upper_MUE      seed
## 1:                     futility          NA        NA        NA 996631745
## 2: stop for futility at interim   0.9965259 0.4075308 0.4330387 996631745
## 3:                     efficacy          NA        NA        NA 579018813
## 4:                         <NA>   0.9987181 0.3409305 0.3518915 579018813
##                                                       file
## 1: sim-2stage_missing_fixC_nonbinding_ar10_power-1_100.rds
## 2: sim-2stage_missing_fixC_nonbinding_ar10_power-1_100.rds
## 3: sim-2stage_missing_fixC_nonbinding_ar10_power-1_100.rds
## 4: sim-2stage_missing_fixC_nonbinding_ar10_power-1_100.rds

library(DelayedGSD)

## * Setting
method <- 3
binding <- FALSE
fixC <- TRUE
ar <- 10
hypo <- "power"
PropForInterim <- 0.5

nGSD <- c(557, 557, 557)

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

## * Example 1
seed1 <- 579018813 

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
thets1 <- df1.sim$t3[ceiling(nGSD*PropForInterim)]

## ** interim
debug1.lmmI <- analyzeData(SelectData(df1.sim,t=thets1[method]),
                           ddf = "nlme", data.decision = sum(df1.sim$t1 <= thets1[method] + 1.50001*14), getinfo = TRUE, trace = TRUE)
debug1.GSDI <- update(debug.plannedB, delta = debug1.lmmI, trace = TRUE)
## summary(debug1.GSDI)

##            GSD with repeated measurements at the interim analysis of stage 1 

## Boundaries and observed statistics 
## stage |         Interim             | Decision       |     Spent          
##       | F-bound E-bound    Stat     |  C-bound Stat  |     alpha      beta
##     1 | 0.35898 2.41613 2.47613 S-E |  1.95996       | 0.0063357 0.0506857
##     2 |                             |                |                    

## Observed and predicted information: 
## stage |  Interim     (%) | Decision     (%) |        n
##     1 | 11.72226 0.50342 | 15.00862 0.64455 | 395.0000
##     2 |                  |                  | 784.6383

## Current ML-estimate of the treatment effect (Z1) 
## estimate      se   lower   upper statistic  df   p.value
##  0.72322 0.29208 0.14948 1.29696   2.47613 541 0.0067932
                                                                   
 ##  * 395 patients in the study: 300 with at least one outcome value 
 ##                               247 with complete data              


## ** decision
debug1.lmmD <- analyzeData(df1.sim[which(df1.sim$t1 <= thets1[method] + 1.50001*14),],
                          ddf = "nlme", getinfo = TRUE, trace = TRUE)
debug1.GSDD <- update(debug1.GSDI, delta = debug1.lmmD, k = 1, type.k = "decision", trace = FALSE)
summary(debug1.GSDD)

## 	    GSD with repeated measurements at the decision analysis of stage 1 

##  Boundaries and observed statistics 
##  stage |         Interim             | Decision           |     Spent          
##        | F-bound E-bound    Stat     |  C-bound    Stat   |     alpha      beta
##      1 | 0.35898 2.41613 2.47613 S-E |  1.95996 1.93293 F | 0.0063357 0.0506857
##      2 |                             |                    |                    

##  Observed and predicted information: 
##  stage |  Interim     (%) | Decision     (%) |   n
##      1 | 11.72226 0.50342 |  16.5819 0.71212 | 395
##      2 |                  |                  |    

##  Current MUE-estimate of the treatment effect (Z1) 
##  estimate   lower   upper  p.value
##   0.34093 0.34093 0.35189 0.998718
                                                                   
##  395 patients in the study: 377 with at least one outcome value 
##                             348 with complete data              

## Advarselsbesked:
## I confint.delayedGSD(x, method = c("ML", "MUE")) :
##   Possibly incorrect evaluation of the MUE lower, estimate.
## Mismatch in the optimization process in term of confidence level: 0.943991718560957, 0.246603974336232. 

## * Example 2
seed2 <- 996631745

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
thets2 <- df2.sim$t3[ceiling(nGSD*PropForInterim)]

## ** interim
debug2.lmmI <- analyzeData(SelectData(df2.sim,t=thets2[method]),
                           ddf = "nlme", data.decision = sum(df2.sim$t1 <= thets2[method] + 1.50001*14), getinfo = TRUE, trace = TRUE)
debug2.GSDI <- update(debug.plannedB, delta = debug2.lmmI, trace = TRUE)
## summary(debug2.GSDI)

##            GSD with repeated measurements at the interim analysis of stage 1 

## Boundaries and observed statistics 
## stage |         Interim             | Decision       |     Spent          
##       | F-bound E-bound    Stat     |  C-bound Stat  |     alpha      beta
##     1 |   0.709 2.24163 0.59188 S-F |  1.95996       | 0.0090704 0.0725629
##     2 |                             |                |                    

## Observed and predicted information: 
## stage |  Interim     (%) | Decision     (%) |        n
##     1 | 14.02574 0.60234 |  17.8922 0.76839 | 404.0000
##     2 |                  |                  | 670.7168

## Current ML-estimate of the treatment effect (Z1) 
## estimate      se    lower   upper statistic  df p.value
##  0.15804 0.26702 -0.36647 0.68255   0.59188 543 0.27709
                                                                   
##  * 404 patients in the study: 308 with at least one outcome value 
##                               241 with complete data              

## ** decision
debug2.lmmD <- analyzeData(df2.sim[which(df2.sim$t1 <= thets2[method] + 1.50001*14),],
                          ddf = "nlme", getinfo = TRUE, trace = TRUE)
debug2.GSDD <- update(debug2.GSDI, delta = debug2.lmmD, k = 1, type.k = "decision", trace = FALSE)
summary(debug2.GSDD)

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
##   0.40752 0.40753 0.43304 0.9965259
                                                                   
##  404 patients in the study: 382 with at least one outcome value 
##                             354 with complete data              

## Advarselsbesked:
## I confint.delayedGSD(x, method = c("ML", "MUE")) :
##   Possibly incorrect evaluation of the MUE lower, estimate.
## Mismatch in the optimization process in term of confidence level: 0.922632457177631, 0.235746748986334. 

## * debug
calcP <- function(delta, estimate, binding){
    FinalPvalue(Info.d=c(17.97582495),
                Info.i=c(14.02573904),
                ck=c(1.95996398),
                ck.unrestricted=c(1.95282645),
                lk=c(0.70900319),
                uk=c(2.24162645),
                reason.interim="futility",
                kMax=2,
                estimate,
                delta=delta,
                method=3,
                bindingFutility=binding,
                cNotBelowFixedc=TRUE,
                continuity.correction=TRUE)
}

## estimate: 0.550969
grid <- data.frame(x=seq(-5,5,length.out=1000))
grid$value.binding <- sapply(grid$x,calcP, estimate = 0, binding = TRUE)
grid$value.nonbinding <- sapply(grid$x,calcP, estimate = 0, binding = FALSE)
plot(grid$x,grid$value.binding)
plot(grid$x,grid$value.nonbinding)
##        0%       25%       50%       75%      100% 
## 0.9855382 1.0000000 1.0000000 1.0000000 1.0000000 

## prob to continue (and therefore have a more extreme result than concluding futility now)
grid$term1.binding <- sapply(grid$x, function(delta){
    mvtnorm::pmvnorm(lower = 0.70900319,  
                     upper = 2.24162645,
                     mean = delta * sqrt(14.02573904),
                     sigma = matrix(1, nrow = 1, ncol = 1))
})
grid$term1.nonbinding <- sapply(grid$x, function(delta){
    mvtnorm::pmvnorm(lower = -Inf,  
                     upper = 2.24162645,
                     mean = delta * sqrt(14.02573904),
                     sigma = matrix(1, nrow = 1, ncol = 1))
})

## prob to stop for futility at interim and conclude a more extreme result 
grid$term2.nonbinding <- grid$term2.binding <- sapply(grid$x, function(delta){
    mvtnorm::pmvnorm(lower = c(-Inf, delta * sqrt(17.97582)),   
                     upper = c(0.70900319,Inf),
                     mean = delta * sqrt(c(14.02573904, 17.97582495)),
                     sigma = matrix(c(1,sqrt(14.02573904/17.97582),sqrt(14.02573904/17.97582),1), nrow = 2, ncol = 2))
})  
      
## prob to stop for efficacy at interim and conclude more extreme result (use min of critval and statistic as stopping for eff and z_k>c_k will result in efficacy conclusion (always more extreme than concluding futility, which is the case when stopping for futility at interim for Method 3 as flips from fut to eff not allowed))
grid$term3.binding <- grid$term3.nonbinding <- sapply(grid$x, function(delta){
    mvtnorm::pmvnorm(lower = c(2.241626, min(1.959964,delta * sqrt(17.97582))),  
                     upper = c(Inf,Inf),
                     mean = delta * sqrt(c(14.02573904, 17.97582495)),
                     sigma = matrix(c(1,sqrt(14.02573904/17.97582),sqrt(14.02573904/17.97582),1), nrow = 2, ncol = 2))
})
range(pmin(grid$term1.nonbinding + grid$term2.nonbinding + grid$term3.nonbinding, 1) - grid$value.nonbinding)


library(ggplot2)
library(reshape2)
gridL <- reshape2::melt(grid, id.vars = "x")
gridL$term <- gsub(".binding|.nonbinding","",as.character(gridL$variable))
gridL$binding <- grepl(".nonbinding",as.character(gridL$variable)) == FALSE

gg.term <- ggplot(gridL, aes(x=x, y = value, color = term, group = variable))
gg.term <- gg.term + geom_line()
gg.term <- gg.term + facet_wrap(~binding, labeller = label_both)
gg.term



FinalPvalue(Info.d=c(17.97582495),
                Info.i=c(14.02573904),
                ck=c(1.95996398),
                ck.unrestricted=c(1.95282645),
                lk=c(0.70900319),
                uk=c(2.24162645),
                reason.interim="futility",
                kMax=2,
                estimate=0,
                delta=5,
                method=2,
                bindingFutility=FALSE,
                cNotBelowFixedc=TRUE,
                continuity.correction=TRUE)
##----------------------------------------------------------------------
### example-debug.R ends here
