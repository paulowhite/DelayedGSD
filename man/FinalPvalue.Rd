% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FinalPvalue.R
\name{FinalPvalue}
\alias{FinalPvalue}
\title{P-value at Decision}
\usage{
FinalPvalue(Info.d, Info.i, ck, lk, uk, sided = 1, kMax, delta = 0, estimate)
}
\arguments{
\item{Info.d}{Information at all decision analyses up to stage where study was stopped (should include the information at the final analysis if stopped at final analysis)}

\item{Info.i}{Information at all interim analyses up to stage where study was stopped}

\item{ck}{decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis)}

\item{lk}{lower bounds up to stage where study was stopped}

\item{uk}{upper bounds up to stage where study was stopped}

\item{sided}{whether the test is 1 or 2-sided}

\item{kMax}{maximum number of analyses}

\item{delta}{true effect under which to calculate the probability (should always be 0 for p-value, only change for calculation of CI)}

\item{estimate}{naive estimate (e.g. using  ML or REML).}
}
\description{
Compute the p-value at the end of the trial.
}
\examples{
library(rpact)

#### example 1 ####
## simulate data for a given design
kMax <- 3 
design <- getDesignGroupSequential(kMax=kMax, sided = 1, alpha = 0.025, beta = 0.2,
                                  informationRates = c(0.5,0.75,1),
                                  typeOfDesign="asKD",
                                  typeBetaSpending="bsKD",gammaA=2,gammaB=2)
efficacyBound <- design$criticalValues
futilityBound <- design$futilityBounds

res <- getDataset(sampleSizes1 = c(40,20),sampleSizes2=c(40,20),means1=c(1.5,3),means2=c(0.5,-1),stDevs1=c(2,1.65),stDevs2=c(2,1.65))

## extract key informations
vec.stage <- res$stages ##  1 1 2 2
vec.mean <- res$overallMeans ##   1.5 0.5 2.0 0.0
vec.std <- res$overallStDevs ##  2.000000 2.000000 2.007307 2.007307
vec.z <- vec.mean/(vec.std/sqrt(vec.n)) ## 4.743416 1.581139 7.717771 0.000000
vec.n <- res$overallSampleSizes ##  40 40 60 60
estimate <- abs(diff(vec.mean[vec.stages==res$.kMax])) ##  2

Info.var <- as.double(1/tapply(vec.std^2/vec.n,vec.stage,sum)) ## 5.0000 7.4455 
statistic <- estimate*sqrt(Info.var[2]) ## 5.457289 

Info.cor <- (rep(1,res$.kMax) \%o\% sqrt(Info.var)) / t(rep(1,res$.kMax) \%o\% sqrt(Info.var))
Info.cor[upper.tri(Info.cor)] <- 0
Info.cor <- Info.cor + t(Info.cor)
diag(Info.cor) <- 1

## p-value by hand
## - more extreme 1: stop at first interim for efficacy
pval1 <- pmvnorm(lower = efficacyBound[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
## - more extreme 1: stop at second interim for efficacy with a larger effect
pval2 <- pmvnorm(lower = c(futilityBound[1], statistic), upper = c(efficacyBound[1],Inf), mean=c(0,0), sigma = Info.cor) 
## global
pval1 + pval2 ## 0.00625 

## p-value using FinalPvalue
FinalPvalue(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 3, estimate = estimate)
FinalPvalue2(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 3, estimate = estimate)

## p-value using rpact

ests <- getAnalysisResults(design, res, normalApproximation=TRUE)

## Final p-value                                : NA, 0.00625, NA 
##  Final CIs (lower)                            : NA, 0.2414, NA 
##  Final CIs (upper)                            : NA, 2.001, NA 
##  Median unbiased estimate                     : NA, 1.121, NA

#### example 2 (see Wassmer 2016, Group sequential and confirmatory adaptive designs in clinical trials, section 4.1, page 87-88) ####
vec.stage <- 1:2 ##  1 1 2 2
vec.std <- c(1,1)
vec.z <- c(1,3)
vec.n <- c(22,44)
estimate <- vec.z/sqrt(vec.n[2])

efficacyBound <- c(2.361,2.361)
futilityBound <- c(-2.361,-2.361)

Info.var <- vec.n
statistic <- estimate*sqrt(Info.var[2]) ## 5.457289 

Info.cor <- (rep(1,length(vec.stage)) \%o\% sqrt(Info.var)) / t(rep(1,length(vec.stage)) \%o\% sqrt(Info.var))
Info.cor[upper.tri(Info.cor)] <- 0
Info.cor <- Info.cor + t(Info.cor)
diag(Info.cor) <- 1

## p-value by hand
## - more extreme 1: stop at first interim for efficacy
pval1 <- pmvnorm(lower = efficacyBound[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
## - more extreme 1: stop at second interim for efficacy with a larger effect
pval2 <- pmvnorm(lower = c(futilityBound[1], statistic[2]), upper = c(efficacyBound[1],Inf), mean=c(0,0), sigma = Info.cor) 
## global
pval1 + pval2 ## 0.009819342 

## p-value using FinalPvalue
FinalPvalue(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate)
FinalPvalue2(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate) ## [1] 0.009819408

## confidence interval using FinalPvalue
FinalPvalue2(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate, delta = 0.0737)
1-FinalPvalue2(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate, delta = 0.729)

FinalCI2(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate) ## [1] 0.07371176 0.72920245

library(microbenchmark)
microbenchmark(optimise = FinalCI2(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate, optimizer = "optimise"),
              uniroot = FinalCI2(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate, optimizer = "uniroot") 
)
## Unit: milliseconds
##      expr      min       lq     mean   median       uq      max neval cld
##  optimise 27.58699 30.60128 32.20530 32.24507 33.47103 44.99762   100   b
##   uniroot 21.37910 23.32777 24.45567 24.29751 25.16585 29.36898   100  a

## p-value by hand
## - more extreme 1: stop at first interim for efficacy
pval1 <- pmvnorm(lower = efficacyBound[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
## - more extreme 1: stop at second interim for efficacy with a larger effect
pval2 <- pmvnorm(lower = c(futilityBound[1], statistic), upper = c(efficacyBound[1],Inf), mean=c(0,0), sigma = Info.cor) 
## global
pval1 + pval2 ## 0.00625 
}
