## * FinalPvalue (documentation)
#' @title P-value at Decision
#' @description  Compute the p-value at the end of the trial.
#' @param Info.d Information at all decision analyses up to stage where study was stopped (should include the information at the final analysis if stopped at final analysis)
#' @param Info.i Information at all interim analyses up to stage where study was stopped
#' @param ck decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis)
#' @param lk lower bounds up to stage where study was stopped
#' @param uk upper bounds up to stage where study was stopped
#' @param kMax maximum number of analyses
#' @param delta true effect under which to calculate the probability (should always be 0 for p-value, only change for calculation of CI)
#' @param estimate naive estimate (e.g. using  ML or REML).
#'
#' @examples
#' library(mvtnorm)
#' 
#' if(require(rpact)){
#'
#' #### example 1 ####
#' ## simulate data for a given design
#' kMax <- 3 
#' design <- getDesignGroupSequential(kMax=kMax, sided = 1, alpha = 0.025, beta = 0.2,
#'                                   informationRates = c(0.5,0.75,1),
#'                                   typeOfDesign="asKD",
#'                                   typeBetaSpending="bsKD",gammaA=2,gammaB=2)
#' efficacyBound <- design$criticalValues
#' futilityBound <- design$futilityBounds
#' 
#' res <- getDataset(sampleSizes1 = c(40,20),sampleSizes2=c(40,20),means1=c(1.5,3),means2=c(0.5,-1),stDevs1=c(2,1.65),stDevs2=c(2,1.65))
#'
#' ## extract key informations
#' vec.stage <- res$stages ##  1 1 2 2
#' vec.mean <- res$overallMeans ##   1.5 0.5 2.0 0.0
#' vec.std <- res$overallStDevs ##  2.000000 2.000000 2.007307 2.007307
#' vec.n <- res$overallSampleSizes ##  40 40 60 60
#' vec.z <- vec.mean/(vec.std/sqrt(vec.n)) ## 4.743416 1.581139 7.717771 0.000000
#' estimate <- abs(diff(vec.mean[vec.stage==res$.kMax])) ##  2
#' 
#' Info.var <- as.double(1/tapply(vec.std^2/vec.n,vec.stage,sum)) ## 5.0000 7.4455 
#' statistic <- estimate*sqrt(Info.var[2]) ## 5.457289 
#' 
#' Info.cor <- (rep(1,res$.kMax) %o% sqrt(Info.var)) / t(rep(1,res$.kMax) %o% sqrt(Info.var))
#' Info.cor[upper.tri(Info.cor)] <- 0
#' Info.cor <- Info.cor + t(Info.cor)
#' diag(Info.cor) <- 1
#' 
#' ## p-value by hand
#' ## - more extreme 1: stop at first interim for efficacy
#' pval1 <- pmvnorm(lower = efficacyBound[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
#' ## - more extreme 1: stop at second interim for efficacy with a larger effect
#' pval2 <- pmvnorm(lower = c(futilityBound[1], statistic), upper = c(efficacyBound[1],Inf), mean=c(0,0), sigma = Info.cor) 
#' ## global
#' pval1 + pval2 ## 0.00625 
#' 
#' ## p-value using FinalPvalue
#' FinalPvalue(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 3, estimate = estimate)
#' 
#' ## p-value using rpact
#' 
#' ests <- getAnalysisResults(design, res, normalApproximation=TRUE)
#'
#' ## Final p-value                                : NA, 0.00625, NA 
#' ##  Final CIs (lower)                            : NA, 0.2414, NA 
#' ##  Final CIs (upper)                            : NA, 2.001, NA 
#' ##  Median unbiased estimate                     : NA, 1.121, NA
#'
#' #### example 2 (see Wassmer 2016, Group sequential and confirmatory adaptive designs in clinical trials, section 4.1, page 87-88) ####
#' vec.stage <- 1:2 ##  1 1 2 2
#' vec.std <- c(1,1)
#' vec.z <- c(1,3)
#' vec.n <- c(22,44)
#' estimate <- vec.z/sqrt(vec.n[2])
#'
#' efficacyBound <- c(2.361,2.361)
#' futilityBound <- c(-2.361,-2.361)
#' 
#' Info.var <- vec.n
#' statistic <- estimate*sqrt(Info.var[2]) ## 5.457289 
#' 
#' Info.cor <- (rep(1,length(vec.stage)) %o% sqrt(Info.var)) / t(rep(1,length(vec.stage)) %o% sqrt(Info.var))
#' Info.cor[upper.tri(Info.cor)] <- 0
#' Info.cor <- Info.cor + t(Info.cor)
#' diag(Info.cor) <- 1
#' 
#' ## p-value by hand
#' ## - more extreme 1: stop at first interim for efficacy
#' pval1 <- pmvnorm(lower = efficacyBound[1], upper = Inf, mean=0, sigma= Info.cor[1,1,drop=FALSE]) 
#' ## - more extreme 1: stop at second interim for efficacy with a larger effect
#' pval2 <- pmvnorm(lower = c(futilityBound[1], statistic[2]), upper = c(efficacyBound[1],Inf), mean=c(0,0), sigma = Info.cor) 
#' ## global
#' pval1 + pval2 ## 0.009819342 
#' 
#' ## p-value using FinalPvalue
#' FinalPvalue(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate)
#' 
#' ## confidence interval using FinalPvalue
#' FinalCI(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate) ## [1] 0.07371176 0.72920245
#'
#' ## library(microbenchmark)
#' ## microbenchmark(optimise = FinalCI(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate, optimizer = "optimise"),
#' ##               uniroot = FinalCI(Info.d = Info.var, Info.i = Info.var, ck = efficacyBound, lk = futilityBound, uk = efficacyBound, kMax = 4, estimate = estimate, optimizer = "uniroot") 
#' ## )
#' ## Unit: milliseconds
#' ##      expr      min       lq     mean   median       uq      max neval cld
#' ##  optimise 27.58699 30.60128 32.20530 32.24507 33.47103 44.99762   100   b
#' ##   uniroot 21.37910 23.32777 24.45567 24.29751 25.16585 29.36898   100  a
#' 
#' }

## * FinalPvalue (code)
#' @export
FinalPvalue <- function(Info.d,  
                        Info.i,  
                        ck,   
                        lk,  
                        uk,  
                        kMax, 
                        delta=0,  
                        estimate){
    
    requireNamespace("mvtnorm")
  
  ## message("the method assumes that positive effects are good")
  k <- length(Info.d)
  I_all <- as.vector(rbind(Info.i, Info.d))  ## vectorize in the following order: |interim 1| |decision 1| |interim 2| ... |interim k| |decision k| |decision final|
  index_tempo <- as.vector(rbind(rep(1,length(Info.i)), rep(2,length(Info.d))))
  index_interim <- which(index_tempo==1) ## 1, 3, ...
  index_decision <- which(index_tempo==2) ## 2, 4, ...
  index_final <- index_decision[kMax] ## NA when the study is stopped before final

  m <- length(I_all)
  sigmaZm <- diag(1,m)
  for(i in 1:m){
    for(j in i:m){
      sigmaZm[i,j] <- sqrt(I_all[i]/I_all[j])
      sigmaZm[j,i] <- sqrt(I_all[i]/I_all[j])
    }
  }
  
  theta <- delta * sqrt(I_all)
  statistic <- estimate * sqrt(Info.d)
  
  index_infoPB <- which(Info.d > Info.d[k])  #find decision information levels that are higher than final information levels.
  
  pval <- 0
  if(statistic[k] >= ck[k]){  #In case efficacy is concluded, only efficacious results at earlier IA are more extreme
      for(i in 1:k){ ## i <- 1

          if(i==1){
              iLk <- NULL
              iUk <- NULL
          }else{
              iLk <- lk[1:(i-1)]
              iUk <- uk[1:(i-1)]
          }

          if(i==k){
              if(i==kMax){ ## is it the final analysis (interim 1:(k-1), final)
                  iIndex <- c(index_interim, index_final)
                  ## prob to continue until final analysis and obtain a higher result than the observed
                  pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,statistic[i]),  
                                                  upper = c(iUk,Inf),
                                                  mean = theta[iIndex],
                                                  sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
              }else{ ## or a decision analysis (interim 1:i, decision i)
                  iIndex <- c(index_interim[1:i], index_decision[i])

                  ## prob to stop for eff at analysis k and conclude more extreme result
                  pval <- pval +  mvtnorm::pmvnorm(lower = c(iLk,uk[i],statistic[i]),  
                                                   upper = c(iUk,Inf,Inf),
                                                   mean = theta[iIndex],
                                                   sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
                  ## prob to stop for fut at analysis k and switch to eff with more extreme result
                  pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,statistic[i]),   
                                                  upper = c(iUk,lk[i],Inf),
                                                  mean = theta[iIndex],
                                                  sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
              }
          } else if(i %in% index_infoPB){
              next #results for larger information levels are never more extreme than smaller information levels with results above c
          } else { ## or an interim analysis where we continued (interim 1:i, decision i)
              iIndex <- c(index_interim[1:i], index_decision[i])

                  ## prob to stop for eff at analysis i and conclude eff
                  pval <- pval +  mvtnorm::pmvnorm(lower = c(iLk,uk[i],ck[i]),  
                                                   upper = c(iUk,Inf,Inf),
                                                   mean = theta[iIndex],
                                                   sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
                  ## prob to stop for fut at analysis i and switch to eff
                  pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,ck[i]),   
                                                  upper = c(iUk,lk[i],Inf),
                                                  mean = theta[iIndex],
                                                  sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
          }
      }
  } else {  #In case futility is concluded, calculate 1-obtaining a less extreme result, i.e. 1-prob of concluding futility at earlier IA
      for(i in 1:k){

          if(i==1){
              iLk <- NULL
              iUk <- NULL
          }else{
              iLk <- lk[1:(i-1)]
              iUk <- uk[1:(i-1)]
          }

          if(i==k){
              if (i==kMax){## is it the final analysis (interim 1:(k-1), final)
                  iIndex <- c(index_interim, index_final)
                  pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf),  #prob to continue until final analysis and obtain a lower result than the observed
                                         upper = c(iUk,statistic[i]),
                                         mean = theta[iIndex],
                                         sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
              }else{## or the decision analysis (interim 1:i, decision i)
                  iIndex <- c(index_interim[1:i], index_decision[i])

                  pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,uk[i],-Inf),  #prob to stop for eff at analysis k and conclude less extreme result
                                                  upper = c(iUk,Inf,statistic[i]),
                                                  mean = theta[iIndex],
                                                  sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
                  pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,-Inf), #prob to stop for fut at analysis k and conclude less extreme result
                                                  upper = c(iUk,lk[i],statistic[i]),
                                                  mean = theta[iIndex],
                                                  sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
              }
          } else if (i %in% index_infoPB){ #results for larger information levels are always more extreme if observed result < c
                                        #i cannot be kMax in this case, since Info.d[k] would be Info.d[kMax] in the definition of the index
              next
          } else {## or an interim analysis where we continued (interim 1:i, decision i)
              iIndex <- c(index_interim[1:i], index_decision[i])

              pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,uk[i],-Inf), #prob to stop for eff at analysis i and switch to fut
                                              upper = c(iUk,Inf,ck[i]),
                                              mean = theta[iIndex],
                                              sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
              pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,-Inf), #prob to stop for fut at analysis i and conclude fut
                                              upper = c(iUk,lk[i],ck[i]),
                                              mean = theta[iIndex],
                                              sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
          }
      }
      pval <- 1-pval
  }
  
  return(pval)
}

