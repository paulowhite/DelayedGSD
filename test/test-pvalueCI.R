### test-pvalueCI.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 27 2022 (09:32) 
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

library(testthat)
library(rpact)
set.seed(10)

context("Computing p-value/CIs/MUestimate after stopping. \n")

## * rpact
test_that("Compare p-value/CIs/MUestimate to rpact in the non-longitudinal case",{
    kMax <- 3 
    design <- rpact::getDesignGroupSequential(kMax=kMax, sided = 1, alpha = 0.025, beta = 0.2,
                                              informationRates = c(0.5,0.75,1),
                                              typeOfDesign="asKD",
                                              typeBetaSpending="bsKD",gammaA=2,gammaB=2)
    efficacyBound <- design$criticalValues
    futilityBound <- design$futilityBounds
 
    res <- rpact::getDataset(sampleSizes1 = c(40,20),sampleSizes2=c(40,20),means1=c(1.5,3),means2=c(0.5,-1),stDevs1=c(2,1.65),stDevs2=c(2,1.65))
    ests <- rpact::getAnalysisResults(design, res, normalApproximation=TRUE)

    ## extract key informations
    vec.stage <- res$stages ##  1 1 2 2
    vec.mean <- res$overallMeans ##   1.5 0.5 2.0 0.0
    vec.std <- res$overallStDevs ##  2.000000 2.000000 2.007307 2.007307
    vec.n <- res$overallSampleSizes ##  40 40 60 60
    vec.z <- vec.mean/(vec.std/sqrt(vec.n)) ## 4.743416 1.581139 7.717771 0.000000
    estimate <- abs(diff(vec.mean[vec.stage==res$.kMax])) ##  2
 
    Info.var <- as.double(1/tapply(vec.std^2/vec.n,vec.stage,sum)) ## 5.0000 7.4455 
    statistic <- estimate*sqrt(Info.var[2]) ## 5.457289 
 
    Info.cor <- (rep(1,res$.kMax) %o% sqrt(Info.var)) / t(rep(1,res$.kMax) %o% sqrt(Info.var))
    Info.cor[upper.tri(Info.cor)] <- 0
    Info.cor <- Info.cor + t(Info.cor)
    diag(Info.cor) <- 1
 
    ## p-value by hand
    ## - more extreme 1: stop at first interim for efficacy
    pval1 <- pmvnorm(lower = efficacyBound[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
    ## - more extreme 1: stop at second interim for efficacy with a larger effect
    pval2 <- pmvnorm(lower = c(futilityBound[1], statistic), upper = c(efficacyBound[1],Inf), mean=c(0,0), sigma = Info.cor) 
    ## global
    expect_equal(as.double(pval1 + pval2), 0.00625, tol = 1e-6)
 
    ## p-value using FinalPvalue
    test.p <- FinalPvalue(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 3, estimate = estimate,
                          futility2efficacy = TRUE, bindingFutility = TRUE)
    expect_equal(as.double(test.p), 0.00625)
    expect_equal(as.double(test.p), as.data.frame(ests)[2,"Final p-value"])
    
    ## confidence interval using FinalPvalue
    test.ci <- FinalCI(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate,
                       bindingFutility = TRUE,futility2efficacy=TRUE)
    expect_equal(as.double(test.ci), c(0.53413706, 1.99336202), tol = 1e-4)
    as.data.frame(ests)[2,"Final CI (lower)"] ## 0.2413634
    as.data.frame(ests)[2,"Final CI (upper)"] ## 2.000625

    ## median unbiased estimate
    test.beta <- FinalEstimate(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate, bindingFutility = TRUE, futility2efficacy=TRUE)
    expect_equal(as.double(test.beta), 1.26709147, tol = 1e-6)
    as.data.frame(ests)[2,"Median unbiased estimate"] ## 1.121087
}) 

## * Wassmer 2016 (Group sequential and confirmatory adaptive designs in clinical trials, section 4.1, page 87-88)

test_that("p-value/CIs/MUestimate with four stage design stopped at interim 2",{
 vec.stage <- 1:2 ##  1 1 2 2
 vec.std <- c(1,1)
 vec.z <- c(1,3)
 vec.n <- c(22,44)
 estimate <- vec.z/sqrt(vec.n[2])

 Info.var <- vec.n
 statistic <- estimate*sqrt(Info.var[2]) ## 5.457289 

 Info.cor <- (rep(1,length(vec.stage)) %o% sqrt(Info.var)) / t(rep(1,length(vec.stage)) %o% sqrt(Info.var))
 Info.cor[upper.tri(Info.cor)] <- 0
 Info.cor <- Info.cor + t(Info.cor)
 diag(Info.cor) <- 1

 ## ** O Brien and Fleming
 ## alpha=0.05, power=0.8, delta = 0.31
 
 uk.Obrien <- c(4.049/sqrt(1:2))
 lk.Obrien <- c(-4.049/sqrt(1:2))
 ck.Obrien <- uk.Obrien
 
 ## p-value by hand
 ## - more extreme 1: stop at first interim for efficacy
 pval1.Obrien <- pmvnorm(lower = uk.Obrien[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
 ## - more extreme 1: stop at second interim for efficacy with a larger effect
 pval2.Obrien <- pmvnorm(lower = c(lk.Obrien[1], statistic[2]), upper = c(uk.Obrien[1],Inf), mean=c(0,0), sigma = Info.cor) 
 ## global
 expect_equal(round(as.double(pval1.Obrien + pval2.Obrien),4), 0.0014)
 
 ## p-value using FinalPvalue
 pval.Obrien <- FinalPvalue(Info.d = Info.var, Info.i = Info.var,
                            ck = ck.Obrien, lk = lk.Obrien, uk = uk.Obrien, kMax = 4, estimate = estimate, bindingFutility = TRUE, futility2efficacy = TRUE)
 expect_equal(as.double(pval.Obrien), 0.001362486, tol = 1e-6)
 
 ## confidence interval using FinalPvalue
 CI.Obrien <- FinalCI(Info.d = Info.var, Info.i = Info.var, ck = ck.Obrien, lk = lk.Obrien, uk = uk.Obrien, kMax = 4, estimate = estimate, bindingFutility = TRUE, futility2efficacy = TRUE)
 expect_equal(round(as.double(CI.Obrien),3), c(0.157, 0.748), tol = 1e-3)

 ## median unbiased estimate
 estimate.Obrien <- FinalEstimate(Info.d = Info.var, Info.i = Info.var, ck = ck.Obrien, lk = lk.Obrien, uk = uk.Obrien, kMax = 4, estimate = estimate, bindingFutility = TRUE, futility2efficacy = TRUE)
 expect_equal(round(as.double(estimate.Obrien),3), 0.452, tol = 1e-3)

 ## ** Pocock case
 uk.Pocock <- rep(2.361,2)
 lk.Pocock <- rep(-2.361,2)
 ck.Pocock <- rep(2.361,2)

 ## p-value by hand
 ## - more extreme 1: stop at first interim for efficacy
 pval1.Pocock <- pmvnorm(lower = uk.Pocock[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
 ## - more extreme 1: stop at second interim for efficacy with a larger effect
 pval2.Pocock <- pmvnorm(lower = c(lk.Pocock[1], statistic[2]), upper = c(uk.Pocock[1],Inf), mean=c(0,0), sigma = Info.cor) 
 ## global
 expect_equal(round(as.double(pval1.Pocock + pval2.Pocock),4), 0.0098)
 
 ## p-value using FinalPvalue
 pval.Pocock <- FinalPvalue(Info.d = Info.var, Info.i = Info.var, ck = ck.Pocock, lk = lk.Pocock, uk = uk.Pocock, kMax = 4, estimate = estimate, bindingFutility = TRUE, futility2efficacy = TRUE)
 expect_equal(as.double(pval.Pocock), 0.009819342, tol = 1e-6) 
 
 ## confidence interval using FinalPvalue
 CI.Pocock <- FinalCI(Info.d = Info.var, Info.i = Info.var, ck = ck.Pocock, lk = lk.Pocock, uk = uk.Pocock, kMax = 4, estimate = estimate, bindingFutility = TRUE, futility2efficacy = TRUE)
 expect_equal(round(as.double(CI.Pocock),3), c(0.074, 0.729), tol = 1e-3)

 ## median unbiased estimate
 estimate.Pocock <- FinalEstimate(Info.d = Info.var, Info.i = Info.var, ck = ck.Pocock, lk = lk.Pocock, uk = uk.Pocock, kMax = 4, estimate = estimate, bindingFutility = TRUE, futility2efficacy = TRUE)
 expect_equal(round(as.double(estimate.Pocock),3), 0.419, tol = 1e-3)
})

## * Paul's binding test
test_that("Check consistency between p-value and boundary",{
    #Method 1
    #for binding futility rule
    myBound <- CalcBoundaries(kMax=3,
                              alpha=0.025, 
                              beta=0.1, 
                              InfoR.i=c(3.5,6.75,12)/12,
                              rho_alpha=1.345,
                              rho_beta=1.345,
                              method=1, ## has been changed from 2 to 1
                              cNotBelowFixedc=FALSE,
                              bindingFutility=TRUE,
                              delta=1,
                              InfoR.d=c(5.5,8.75)/12)

            
    test <- FinalPvalue(Info.d = myBound$Info.d,
                        Info.i = myBound$Info.i,
                        ck = myBound$ck,
                        lk = myBound$lk,
                        uk = myBound$uk,
                        kMax = myBound$kMax,
                        estimate = myBound$uk[myBound$kMax]/sqrt(myBound$Info.max),
                        futility2efficacy = TRUE,
                        bindingFutility = TRUE)

    expect_equal(as.double(test), 0.025, tol = 1e-3)

    #for non-binding futility rule
    myBound <- CalcBoundaries(kMax=3,
                              alpha=0.025, 
                              beta=0.1, 
                              InfoR.i=c(3.5,6.75,12)/12,
                              rho_alpha=1.345,
                              rho_beta=1.345,
                              method=1, 
                              cNotBelowFixedc=FALSE,
                              bindingFutility=FALSE,
                              delta=1,
                              InfoR.d=c(5.5,8.75)/12)

            
    test <- FinalPvalue(Info.d = myBound$Info.d,
                        Info.i = myBound$Info.i,
                        ck = myBound$ck,
                        lk = myBound$lk,
                        uk = myBound$uk,
                        kMax = myBound$kMax,
                        estimate = myBound$uk[myBound$kMax]/sqrt(myBound$Info.max),
                        futility2efficacy = TRUE,
                        bindingFutility = FALSE)

    expect_equal(as.double(test), 0.025, tol = 1e-3)
    
    #Method 2
    #for binding futility rule
    myBound <- CalcBoundaries(kMax=3,
                              alpha=0.025, 
                              beta=0.1, 
                              InfoR.i=c(3.5,6.75,12)/12,
                              rho_alpha=1.345,
                              rho_beta=1.345,
                              method=2, ## has been changed from 2 to 1
                              cNotBelowFixedc=FALSE,
                              bindingFutility=TRUE,
                              delta=1,
                              InfoR.d=c(5.5,8.75)/12)
    
    
    test <- FinalPvalue(Info.d = myBound$Info.d,
                        Info.i = myBound$Info.i,
                        ck = myBound$ck,
                        lk = myBound$lk,
                        uk = myBound$uk,
                        kMax = myBound$kMax,
                        estimate = myBound$uk[myBound$kMax]/sqrt(myBound$Info.max),
                        futility2efficacy = TRUE,
                        bindingFutility = TRUE)
    
    expect_equal(as.double(test), 0.025, tol = 1e-3)
    
    #for non-binding futility rule
    myBound <- CalcBoundaries(kMax=3,
                              alpha=0.025, 
                              beta=0.1, 
                              InfoR.i=c(3.5,6.75,12)/12,
                              rho_alpha=1.345,
                              rho_beta=1.345,
                              method=2, 
                              cNotBelowFixedc=FALSE,
                              bindingFutility=FALSE,
                              delta=1,
                              InfoR.d=c(5.5,8.75)/12)
    
    
    test <- FinalPvalue(Info.d = myBound$Info.d,
                        Info.i = myBound$Info.i,
                        ck = myBound$ck,
                        lk = myBound$lk,
                        uk = myBound$uk,
                        kMax = myBound$kMax,
                        estimate = myBound$uk[myBound$kMax]/sqrt(myBound$Info.max),
                        futility2efficacy = TRUE,
                        bindingFutility = FALSE)
    
    expect_equal(as.double(test), 0.025, tol = 1e-3)
    
    #Method 3
    #for binding futility rule
    myBound <- CalcBoundaries(kMax=3,
                              alpha=0.025, 
                              beta=0.1, 
                              InfoR.i=c(3.5,6.75,12)/12,
                              rho_alpha=1.345,
                              rho_beta=1.345,
                              method=3, 
                              cNotBelowFixedc=TRUE,
                              bindingFutility=TRUE,
                              delta=1,
                              InfoR.d=c(5.5,8.75)/12)
    
    
    test <- FinalPvalue(Info.d = myBound$Info.d,
                        Info.i = myBound$Info.i,
                        ck = myBound$ck,
                        lk = myBound$lk,
                        uk = myBound$uk,
                        kMax = myBound$kMax,
                        estimate = myBound$uk[myBound$kMax]/sqrt(myBound$Info.max),
                        futility2efficacy = FALSE,
                        bindingFutility = TRUE)
    
    expect_equal(as.double(test), 0.025, tol = 1e-3)
    
    #for non-binding futility rule
    myBound <- CalcBoundaries(kMax=3,
                              alpha=0.025, 
                              beta=0.1, 
                              InfoR.i=c(3.5,6.75,12)/12,
                              rho_alpha=1.345,
                              rho_beta=1.345,
                              method=3, 
                              cNotBelowFixedc=TRUE,
                              bindingFutility=FALSE,
                              delta=1,
                              InfoR.d=c(5.5,8.75)/12)
    
    
    test <- FinalPvalue(Info.d = myBound$Info.d,
                        Info.i = myBound$Info.i,
                        ck = myBound$ck,
                        lk = myBound$lk,
                        uk = myBound$uk,
                        kMax = myBound$kMax,
                        estimate = myBound$uk[myBound$kMax]/sqrt(myBound$Info.max),
                        futility2efficacy = FALSE,
                        bindingFutility = FALSE)
    
    expect_equal(as.double(test), 0.025, tol = 1e-3)

})
##----------------------------------------------------------------------
### test-pvalueCI.R ends here

