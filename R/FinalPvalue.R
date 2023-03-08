## * FinalPvalue (documentation)
#' @title P-value at Decision
#' @description  Compute the p-value at the end of the trial.
#' @param Info.d Information at all decision analyses up to stage where study was stopped (excepted the final analysis)
#' @param Info.i Information at all interim analyses up to stage where study was stopped (should include the information at the final analysis if stopped at final analysis)
#' @param ck,ck.unrestricted decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis).
#' ck is possibly with restriction (when cNotBelowFixedc=TRUE) and ck.unrestricted always without.
#' @param lk lower bounds up to stage where study was stopped (should include final boundary if stopped at final analysis)
#' @param uk upper bounds up to stage where study was stopped (should include final boundary if stopped at final analysis)
#' @param kMax maximum number of analyses
#' @param delta true effect under which to calculate the probability (should always be 0 for p-value, only change for calculation of CI)
#' @param estimate naive estimate (e.g. using  ML or REML).
#' @param method  method 1, 2 or 3
#' @param bindingFutility [logical] whether the futility stopping rule is binding.
#' @param cNotBelowFixedc [logical] whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
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
#' 
#' 
#' ##Example to show that the p-value will be 0.025 if the final test statistic is exactly on the boundary (requires cNotBelowFixedc=F and in case of non-binding futility that the info at decision is same as at interim)
#' myBound <- CalcBoundaries(kMax=3,
#'                           sided=1,
#'                           alpha=0.025,  
#'                           beta=0.2,  
#'                           InfoR.i=c(1/3,2/3,1),
#'                           rho_alpha=2,
#'                           rho_beta=2,
#'                           method=1,       #works both for method=1 and method=2
#'                           cNotBelowFixedc=FALSE,
#'                           bindingFutility=FALSE,
#'                           delta=1,
#'                           InfoR.d=c(1/3,2/3))
#' 
#' FinalPvalue(myBound$Info.d,
#'             myBound$Info.i,
#'             myBound$ck,
#'             myBound$lk,
#'             myBound$uk,
#'             myBound$kMax,
#'             estimate = myBound$uk[myBound$kMax]/sqrt(myBound$Info.max),
#'             bindingFutility=F,
#'             futility2efficacy=T) 

## * FinalPvalue (code)
#' @export
FinalPvalue <- function(Info.d,  
                        Info.i,  
                        ck,
                        ck.unrestricted,   
                        lk,  
                        uk,  
                        kMax, 
                        delta=0,  
                        estimate,
                        method,  
                        bindingFutility,
                        cNotBelowFixedc){
    
  
  
    requireNamespace("mvtnorm")
  
    ## message("the method assumes that positive effects are good")
    k <- length(Info.i)
    if(k==kMax){
        I_all <- c(as.vector(rbind(Info.i[-kMax], Info.d)),Info.i[kMax])  ## vectorize in the following order: |interim 1| |decision 1| |interim 2| ... |interim k| |decision k| |decision final|
        index_tempo <- as.vector(rbind(rep(1,length(Info.i[-kMax])), rep(2,length(Info.d))))
        index_interim <- which(index_tempo==1) ## 1, 3, ...
        index_decision <- which(index_tempo==2) ## 2, 4, ...
        index_final <- length(I_all)
        critval <- uk[kMax]
    }else{
        I_all <- as.vector(rbind(Info.i, Info.d))  ## vectorize in the following order: |interim 1| |decision 1| |interim 2| ... |interim k| |decision k|
        index_tempo <- as.vector(rbind(rep(1,length(Info.i)), rep(2,length(Info.d))))
        index_interim <- which(index_tempo==1) ## 1, 3, ...
        index_decision <- which(index_tempo==2) ## 2, 4, ...
        index_final <- NA
        critval <- ck
    }

    m <- length(I_all)
    sigmaZm <- diag(1,m)
    for(i in 1:m){
        for(j in i:m){
            sigmaZm[i,j] <- sqrt(I_all[i]/I_all[j])
            sigmaZm[j,i] <- sqrt(I_all[i]/I_all[j])
        }
    }
  
    theta <- delta * sqrt(I_all)
    if(k==kMax){
        statistic <- estimate * sqrt(Info.i[kMax])
        index_infoPB <- which(Info.d > Info.d[k-1])  #find decision information levels that are higher than final information levels.
        index_infoDecr <- which(Info.d < Info.i[1:(kMax-1)]) #find decision information levels that are lower than corresponding interim information levels
        Info.d[index_infoDecr] <- Info.i[index_infoDecr]+Info.i[index_infoDecr]/10000  #set decision analysis information levels to interim information levels in case of decreasing information
    }else{
        statistic <- estimate * sqrt(Info.d[k])
        if(cNotBelowFixedc && (statistic < ck[k]) && (statistic >= ck.unrestricted[k])){
            statistic <- ck.unrestricted[k]-1e-5 ## correction to ensure consistency between rejection and p-value (we do not reject in this interval)
        }
        index_infoPB <- which(Info.d > Info.d[k])  #find decision information levels that are higher than final information levels.
        index_infoDecr <- which(Info.d < Info.i) #find decision information levels that are lower than corresponding interim information levels
        Info.d[index_infoDecr] <- Info.i[index_infoDecr]+Info.i[index_infoDecr]/10000 #set decision analysis information levels to interim information levels in case of decreasing information
    }
    lk_orig <- lk  #need to store original bounds, even in non-binding case to account for probability of switching from futility to efficacy at decision analysis k
    
    if(!bindingFutility){
        lk[1:min(k,kMax-1)] <- -Inf
    }
    pval <- 0
    if(statistic >= critval){  #In case efficacy is concluded, only efficacious results at earlier IA are more extreme
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
                    pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,statistic),  
                                                    upper = c(iUk,Inf),
                                                    mean = theta[iIndex],
                                                    sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
                }else{ ## or a decision analysis (interim 1:i, decision i)
                  iIndex <- c(index_interim[1:i], index_decision[i])  
                  
                  if(method%in%c(1,2)){

                      ## prob to stop for eff at analysis k and conclude more extreme result
                      pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,uk[i],statistic),  
                                                       upper = c(iUk,Inf,Inf),
                                                       mean = theta[iIndex],
                                                       sigma = sigmaZm[iIndex,iIndex,drop=FALSE])

                      ## prob to stop for fut at analysis k and switch to eff with more extreme result
                      pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,statistic),   
                                                      upper = c(iUk,lk_orig[i],Inf),
                                                      mean = theta[iIndex],
                                                      sigma = sigmaZm[iIndex,iIndex,drop=FALSE])
                    } else if (method%in%3){
                      ## prob to stop for eff at analysis k and conclude more extreme result
                      pval <- pval +  mvtnorm::pmvnorm(lower = c(iLk,uk[i],statistic),  
                                                       upper = c(iUk,Inf,Inf),
                                                       mean = theta[iIndex],
                                                       sigma = sigmaZm[iIndex,iIndex,drop=FALSE])
                    }
                    
                 
              }
          } else if(i %in% index_infoPB){
              next #results for larger information levels are never more extreme than smaller information levels with results above c
          } else { ## or an interim analysis where we continued (interim 1:i, decision i)
            
            if(i %in% index_infoDecr & delta==0){ #in case of decreasing information from interim to decision, for calculation of p-value (under H0) only
            #we will ignore the decision analysis and simply calculate the probability of stopping for efficacy. This is correct for Methods 1 and 2 if c can be chosen freely. It is conservative otherwise as it is more likely to flip from efficacy to futility, hence the p-value will be larger than it should be
                
                  iIndex <- c(index_interim[1:i])
                  
                  ## prob to stop for eff at analysis i
                  pval <- pval +  mvtnorm::pmvnorm(lower = c(iLk,uk[i]),  
                                                   upper = c(iUk,Inf),
                                                   mean = theta[iIndex],
                                                   sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
  
            } else {
              if(method%in%c(1,2)){
                if(delta!=0){  #if c has to be at least the critval of a fixed trial, or if we are using this function under delta!=0, then the switching probabilities are not equal
                  iIndex <- c(index_interim[1:i], index_decision[i])
                  
                  ## prob to stop for eff at analysis i and conclude eff
                  pval <- pval +  mvtnorm::pmvnorm(lower = c(iLk,uk[i],ck.unrestricted[i]),  
                                                   upper = c(iUk,Inf,Inf),
                                                   mean = theta[iIndex],
                                                   sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
                  
                  ## prob to stop for fut at analysis i and switch to eff
                  pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,ck.unrestricted[i]),   
                                                  upper = c(iUk,lk_orig[i],Inf),
                                                  mean = theta[iIndex],
                                                  sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
                  
                } else { #if c can be below the critval of a fixed trial, make use of knowledge that switching probs are equal (automatically solves issue with decreasing info at decision as this info is not needed)
                  iIndex <- c(index_interim[1:i])
                  
                  ## prob to stop for eff at analysis i
                  pval <- pval +  mvtnorm::pmvnorm(lower = c(iLk,uk[i]),  
                                                   upper = c(iUk,Inf),
                                                   mean = theta[iIndex],
                                                   sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
                  
                }
              } else if(method%in%3){
                iIndex <- c(index_interim[1:i], index_decision[i])
                ## prob to stop for eff at analysis i and conclude eff
                pval <- pval +  mvtnorm::pmvnorm(lower = c(iLk,uk[i],ck.unrestricted[i]),  
                                                 upper = c(iUk,Inf,Inf),
                                                 mean = theta[iIndex],
                                                 sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
              }
            }
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
                                                    upper = c(iUk,statistic),
                                                    mean = theta[iIndex],
                                                    sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
                }else{## or the decision analysis (interim 1:i, decision i)
                    iIndex <- c(index_interim[1:i], index_decision[i])
                    pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,uk[i],-Inf),  #prob to stop for eff at analysis k and conclude less extreme result
                                                    upper = c(iUk,Inf,statistic),
                                                    mean = theta[iIndex],
                                                    sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
                    pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,-Inf), #prob to stop for fut at analysis k and conclude less extreme result
                                                    upper = c(iUk,lk_orig[i],statistic),
                                                    mean = theta[iIndex],
                                                    sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
                }
            } else if (i %in% index_infoPB){ #results for larger information levels are always more extreme if observed result < c
                                        #i cannot be kMax in this case, since Info.d[k] would be Info.d[kMax] in the definition of the index
                next
            } else {## or an interim analysis where we continued (interim 1:i, decision i)
              
            if(i %in% index_infoDecr & delta==0 & method%in%c(1,2) & bindingFutility){ #in case of decreasing information between interim and decision and we are calculating a p-value (under H0). Only in case of binding futility, otherwise we assume that we ignore option to stop for futility
              iIndex <- c(index_interim[1:i])
              
              pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf), #prob to stop for fut at analysis i. This equals the probability to conclude futility as the switching probabilities are equal
                                              upper = c(iUk,lk_orig[i]),
                                              mean = theta[iIndex],
                                              sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
            } else {
              iIndex <- c(index_interim[1:i], index_decision[i])
              
              pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,uk[i],-Inf), #prob to stop for eff at analysis i and switch to fut
                                              upper = c(iUk,Inf,ck.unrestricted[i]),
                                              mean = theta[iIndex],
                                              sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
              if(method%in%c(1,2) & bindingFutility){  #in case of non-binding rule we will ignore the option to stop for futility, so this part should not be added
                pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,-Inf), #prob to stop for fut at analysis i and conclude fut
                                                upper = c(iUk,lk_orig[i],ck.unrestricted[i]),
                                                mean = theta[iIndex],
                                                sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
              } else if(method%in%3 & bindingFutility){
                pval <- pval + mvtnorm::pmvnorm(lower = c(iLk,-Inf,-Inf), #prob to stop for fut at analysis i and conclude fut (same as stopping for futility for method 3)
                                                upper = c(iUk,lk_orig[i],Inf),
                                                mean = theta[iIndex],
                                                sigma= sigmaZm[iIndex,iIndex,drop=FALSE])
              }
            }
            
            
              
            }
        }
        pval <- 1-pval
    }
  
  return(pval)
}

