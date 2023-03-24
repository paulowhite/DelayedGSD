### example-mismatch-pVal-boundary2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 20 2023 (17:36) 
## Version: 
## Last-Updated: mar 23 2023 (10:53) 
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

library(DelayedGSD)
source("FCT.R")

## * example with crazy p-value
ar <- 5
binding <- FALSE
fixC <- TRUE
delta.factor <- switch("power",
                       "power" = 0.6,
                       "typeI" = 0)

method <- 1
seed <- 187033903

thets <- c(247,249,247)
nGSD <- c(557,561,557)


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
inflationFactor <- debug.plannedB$planned$InflationFactor

df.sim <- GenData(n = max(nGSD), 
                  N.fw = 2,
                  rand.block = c(1,1,0,0),
                  allsd = c(2.5,2.1,2.4),
                  mean0 = c(10,0,0),
                  delta = c(0,0.5,1)*delta.factor,
                  ar = (0.86*2)*2*ar,
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
debug.GSD <- update(debug.plannedB, delta = debug.lmmI, trace = TRUE)

## ** decision
debug.lmmD <- analyzeData(df.sim[which(df.sim$t1 <= thets[method] + 1.50001*14),],
                          ddf = "nlme", getinfo = TRUE, trace = TRUE)
debug.GSD <- update(debug.GSD, delta = debug.lmmD, trace = FALSE)
summary(debug.GSD)

## * example where concluding for futility but with significant p-value
ar <- 10
binding <- FALSE
fixC <- FALSE
delta.factor <- switch("power",
                       "power" = 0.6,
                       "typeI" = 0)

method <- 1
seed <- 520254210
n <- ceiling(2*2*((2.310865/0.6)^2)*(qnorm(1-0.2)-qnorm(0.025))^2)/(1-(0.04807692+0.05769231))

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
inflationFactor <- debug.plannedB$planned$InflationFactor
nGSD <- ceiling(n*inflationFactor)

df.sim <- GenData(n = nGSD, 
                  N.fw = 2,
                  rand.block = c(1,1,0,0),
                  allsd = c(2.5,2.1,2.4),
                  mean0 = c(10,0,0),
                  delta = c(0,0.5,1)*delta.factor,
                  ar = (0.86*2)*2*ar,
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
thets <- c(218,218,218)
nGSD <- c(486,486,486)

## ** interim
SelectData(df.sim,t=thets[method])

debug.lmmI <- analyzeData(SelectData(df.sim,t=thets[method]),
                          ddf = "nlme", data.decision = sum(df.sim$t1 <= thets[method] + 1.50001*14), getinfo = TRUE, trace = TRUE)
debug.GSD <- update(debug.plannedB, delta = debug.lmmI, trace = TRUE)

## ** decision
debug.lmmD <- analyzeData(df.sim[which(df.sim$t1 <= thets[method] + 1.50001*14),],
                          ddf = "nlme", getinfo = TRUE, trace = TRUE)
if(debug.GSD$conclusion["interim",1]=="stop"){
    debug.GSD <- update(debug.GSD, delta = debug.lmmD, trace = FALSE)
}else{
    debug.GSD <- update(debug.GSD, delta = debug.lmmD, k = 1, type.k = "decision", trace = FALSE)
}
debug.GSD$ck
debug.GSD$ck.unrestricted

## ** final
if(debug.GSD$conclusion["interim",1]=="continue"){
    debug.lmmF <- analyzeData(df.sim[1:nGSD[method],],
                              ddf = "nlme", getinfo = TRUE, trace = TRUE)
    debugF.GSD <- update(debug.GSD, delta = debug.lmmF, trace = FALSE)
    summary(debugF.GSD)
    ## confint(debug.GSD)

}    

## * how to get there
if(FALSE){
    library(data.table)
    res2stage <- readRDS("c:/Users/hpl802/Documents/Github/DelayedGSD/Simulations/brice-biostat-cluster/Results-built/res2stage.rds")
    xxx <- res2stage[decision=="futility" & fixC==FALSE & p.value_MUE<0.025,
                     .(scenario,missing,binding,fixC,ar,hypo,method,stage,type,statistic,p.value_ML,p.value_MUE,
                       ck,decision,reason,seed)]
    xxx[,.N,by=c("seed","scenario")][which.max(N)]
    res2stage[seed == 520254210 & scenario == 13 & method == 1, .SD, .SDcols = names(xxx)]
    ##    scenario missing binding  fixC ar  hypo method stage     type statistic p.value_ML p.value_MUE       ck decision              reason      seed
    ## 1:       13    TRUE   FALSE FALSE 10 power      1     1  interim  1.313436         NA          NA       NA continue no boundary crossed 520254210
    ## 2:       13    TRUE   FALSE FALSE 10 power      1     1 decision        NA         NA          NA       NA     <NA>                <NA> 520254210
    ## 3:       13    TRUE   FALSE FALSE 10 power      1     2    final  2.019077 0.02187061  0.02472083 2.038263 futility                <NA> 520254210


    res2stage[seed == 490137693 & scenario == 5, .SD, .SDcols = names(xxx)]
    ##    scenario missing binding fixC ar  hypo method stage     type statistic p.value_ML p.value_MUE       ck decision   reason      seed
    ## 1:        5    TRUE    TRUE TRUE 10 power      1     1  interim  2.542933         NA          NA       NA     stop efficacy 490137693
    ## 2:        5    TRUE    TRUE TRUE 10 power      1     1 decision  1.918589 0.02771629 0.005601226 1.959964 futility     <NA> 490137693
    ## 3:        5    TRUE    TRUE TRUE 10 power      1     2    final        NA         NA          NA       NA     <NA>     <NA> 490137693
    ## 4:        5    TRUE    TRUE TRUE 10 power      2     1  interim  2.542933         NA          NA       NA     stop efficacy 490137693
    ## 5:        5    TRUE    TRUE TRUE 10 power      2     1 decision  1.918589 0.02771629 0.005502683 1.959964 futility     <NA> 490137693
    ## 6:        5    TRUE    TRUE TRUE 10 power      2     2    final        NA         NA          NA       NA     <NA>     <NA> 490137693
    ## 7:        5    TRUE    TRUE TRUE 10 power      3     1  interim  2.542933         NA          NA       NA     stop efficacy 490137693
    ## 8:        5    TRUE    TRUE TRUE 10 power      3     1 decision  1.918589 0.02771629 0.006722260 1.959964 futility     <NA> 490137693
    ## 9:        5    TRUE    TRUE TRUE 10 power      3     2    final        NA         NA          NA       NA     <NA>     <NA> 490137693
}
##----------------------------------------------------------------------
### example-mismatch-pVal-boundary2.R ends here
