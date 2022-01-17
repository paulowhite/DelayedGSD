## * Method2 (code)
updateMethod2 <- function(rho_alpha=2,          # rho parameter of the rho-family spending functions (Kim-DeMets) for alpha
                          rho_beta=2,           # rho parameter of the rho-family spending functions (Kim-DeMets) for beta
                          alpha=0.025, alphaSpent = NULL,         # Type-I error (overall) or spent up to the interim analysis (i.e. cumulative)
                          beta=0.2, betaSpent = NULL,            # Type-II error (overall) or spent at to the interim analysis (i.e. cumulative)
                          Kmax,                 # number of planned analyses (including the final analysis)
                          Info.max=NULL,        # Info.max, i.e. maximum information needed for given beta (type II-error), delta (expected difference), alpha (type I-error) and Kmax. It can be given if it is known. Otherwise it is computed from the  values given for alpha, beta, delta and Kmax.
                          InfoR.i=NULL,         # Expected or observed (wherever possible) information rates at the interim and final analyses 1:Kmax
                          InfoR.d=NULL,         # Expected or observed information rates at all potential decision analyses 1:(Kmax-1)
                          uk = NULL,            # efficacy boundary from the previous interim analyses or planning
                          lk = NULL,            # futility boundary from the previous interim analyses or planning
                          k = NULL, type.k = NULL, ImaxAnticipated = FALSE, # current stage, type of analysis, and conclusion for all previous analyses
                          delta=0,              # expected effect under the alternative (should be on the scale of the test statistc for which If and Info.max relate to one over the variance, e.g. delta=expected log(Hazard ratio))
                          abseps = 1e-06,       # tolerance for precision when finding roots or computing integrals
                          alternative="less",   # greater is for Ho= theta > 0, "less" is for Ho= theta < 0 (note that in Jennison and Turnbull's book chapter (2013) they consider less)
                          binding=TRUE,         # whether the futility boundary is binding
                          Trace=FALSE,          # Used only if Info.max=NULL. Whether to print informations to follow the progression of the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax).
                          nWhileMax=30,         # Used only if Info.max=NULL. Maximum number of steps in the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax)
                          toldiff= 1e-05,       # Used only if Info.max=NULL. Maximum tolerated difference between lower and upper bounds at anaylis Kmax (which souhld be zero), in the root finding algorithm, to find the value of Info.max
                          tolcoef= 1e-04,       # Used only if Info.max=NULL. Maximum tolerated difference before stopping the search (in the root finding algorithm), between two successive values for the multiplier coeficient 'coef' such that Info.max=coef*If  (some values for coef are given in Table 7.6 page 164 Jennison's book. The value of "If" (which stands for Information for Fixed design) corresponds to the information we need if Kmax=1)
                          mycoefMax= 1.2,       # Used only if Info.max=NULL. Upper limit of the interval of values in which we search for the multiplier coeficient 'coef' such that Info.max=coef*If (in the root finding algorithm).
                          mycoefL=1,            # Used only if Info.max=NULL. Lower limit of the interval (see mycoefMax)
                          myseed=2902,           # seed for producing reproducible results. Because we call functions which are based on Monte-Carlo compuation (pmvnorm)
                          cMin=-Inf           # minimun possible value c for the decision analysis, typically that for a fixed sample test (H & J page 10)
                          ){
    ## {{{ set seed


    old <- .Random.seed # to save the current seed
    on.exit( .Random.seed <<- old)
    set.seed(myseed)
    ## }}}
    ## {{{ preliminaries
    ## *** initialize boundaries
    ## lk <- rep(-Inf,Kmax)
    ## uk <- rep(Inf,Kmax)    
    ## initialize

    ## ## *** identify interim analyses that have been skipped due to decreasing information
    ##         index.DI <- which(object$conclusion["reason.interim",]=="decreasing information")
    ##         index.nDI <- setdiff(1:kMax,index.DI)
    ##         if(length(index.DI)>0){
    ##             ## When the interim analysis was skiped because of decreasing information
    ##             ## remove interim and decision
    ##             Info.i <- Info.i[index.nDI]
    ##             Info.d <- Info.d[setdiff(index.nDI,kMax)]
    ##             kMax <- kMax - length(index.DI) ## need to be after the update of InfoR.d as the previous line dependens on kMax
    ##         }



    ##         kMax <- kMax + length(index.DI)

    
    ## *** compute variance-covariance matrix of vector (Z_1,...,Z_k)
    sigmaZk <- diag(1,Kmax)
    for(i in 1:Kmax){
        for(j in i:Kmax){
            sigmaZk[i,j] <- sqrt(InfoR.i[i]/InfoR.i[j])
            sigmaZk[j,i] <- sqrt(InfoR.i[i]/InfoR.i[j])
        }
    }
                                        # compute If (see Jennison book page 87
    ## }}}
    ## {{{ compute vector of means under the alternative H1
  
  #compute information at each analysis
  Info.i <- InfoR.i*Info.max
  Info.d <- InfoR.d*Info.max
  
                                        # compute the mean of the multivariate normal distribution under the alternative H1
    thetheta <- delta*sqrt(Info.i)
    ## }}} 
    ## ** case k=1
    if(k==1){
        ## {{{ case k=1
        if(type.k=="interim"){
            alphaSpent[1] <- ErrorSpend(I=Info.i[1],rho=rho_alpha,beta_or_alpha=alpha,Info.max=Info.max)
            betaSpent[1] <-  ErrorSpend(I=Info.i[1],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)
  
            uk[1] <- qnorm(p=1-alphaSpent[1],mean=0,sd=1)         # compute under the null (Ho)
                                        #thea[1] <- qnorm(p=IncBeta[1],mean=thetheta[1],sd=1)  # compute under the alternative (H1)

                                        #------------
            find.lk <- function(x){
                                        #calculate c corresponding to lk
                ck <- Method1(uk=uk[1],
                              lk=x,
                              Info.i=Info.i[1],
                              Info.d=Info.d[1],
                              Info.max=Info.max,
                              cMin=cMin,
                              ImaxAnticipated=FALSE, ## Should not be here at interim with I>Imax because this case is handled in the parent function (updateBoundary)
                              rho_alpha=rho_alpha,
                              alpha=alpha,
                              bindingFutility = binding)
    
                ## information matrix for first interim and decision analysis
                sigmaZk2 <- matrix(NA,ncol=2,nrow=2)
                sigmaZk2[1,1] <- sigmaZk[1,1]
                sigmaZk2[2,2] <- 1
                sigmaZk2[1,2] <- sigmaZk2[2,1] <- sqrt(Info.i[1]/Info.d[1])
    
                ## probability to conclude futility
                return(pmvnorm(lower = c(uk[1],-Inf),
                               upper = c(Inf,ck),
                               mean=c(thetheta[1],delta*sqrt(Info.d[1])),
                               sigma= sigmaZk2,
                               abseps = abseps) +
                       pmvnorm(lower = c(-Inf,-Inf),
                               upper = c(x,ck),
                               mean=c(thetheta[1],delta*sqrt(Info.d[1])),
                               sigma= sigmaZk2,
                               abseps = abseps) - betaSpent[1])
      
            }
            lk[1] <- uniroot(find.lk,lower=uk[1]-10,upper=uk[1])$root  #dirty solution to use -10 for lower bound
        }
        
        ck <- Method1(uk=uk[1],
                      lk=lk[1],
                      Info.i=Info.i[1],
                      Info.d=Info.d[1],
                      Info.max=Info.max,
                      cMin=cMin,
                      ImaxAnticipated=ImaxAnticipated,
                      rho_alpha=rho_alpha,  
                      alpha=alpha,
                      bindingFutility = binding)
                                        #------------
    }else{
        
        alphaSpentInc <- diff(c(0,alphaSpent))
        betaSpentInc <- diff(c(0,betaSpent))
        
        if(binding){
          TheLowerValues <- lk[1:(k-1)]
        }else{
          TheLowerValues <- rep(-Inf,k-1)
        }
        
        if(type.k=="interim"){

            ## ** Estimate uk
            alphaSpent[k] <- ErrorSpend(I=Info.i[k],rho=rho_alpha,beta_or_alpha=alpha,Info.max=Info.max) 
            alphaSpentInc[k] <- alphaSpent[k] - alphaSpent[(k-1)]   
            betaSpent[k] <- ErrorSpend(I=Info.i[k],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)  
            betaSpentInc[k] <- betaSpent[k] - betaSpent[(k-1)]   
            
            ## {{{ 
            ## {{{ u_k by solving what follows
            uk[k] <- uniroot(function(x){pmvnorm(lower = c(TheLowerValues,x),
                                                 upper = c(uk[1:(k-1)],Inf),
                                                 mean=rep(0,k),
                                                 sigma= sigmaZk[1:k,1:k],
                                                 abseps = abseps) - alphaSpentInc[k]},
                             lower = lk[utils::tail(intersect(which(!is.infinite(lk)),1:(k-1)),1)],  ## last boundary among the k-1 already computed that is not infinite  
                             upper = uk[utils::tail(intersect(which(!is.infinite(uk)),1:(k-1)),1)],
                             ## lower = lk[k-1],
                             ## upper = uk[k-1],
                             tol = abseps)$root

            ## ** Estimate lk
            find.lkk <- function(x){
                ## calculate c corresponding to lk
                ck <- Method1(uk=uk[1:k],
                              lk=c(lk[1:(k-1)],x),
                              Info.i=Info.i[1:k],
                              Info.d=Info.d[k],
                              Info.max=Info.max,
                              cMin=cMin,
                              ImaxAnticipated=FALSE, ## Should not be here at interim with I>Imax because this case is handled in the parent function (updateBoundary)
                              rho_alpha=rho_alpha,
                              alpha=alpha,
                              bindingFutility = binding)
                
                ## information matrix for first interim and decision analysis
                sigmaZk2 <- matrix(NA,ncol=k+1,nrow=k+1)
                sigmaZk2[1:k,1:k] <- sigmaZk[1:k,1:k]
                sigmaZk2[k+1,k+1] <- 1
                sigmaZk2[1:k,k+1] <- sigmaZk2[k+1,1:k] <- sqrt(Info.i[1:k]/Info.d[k])
                
                ## probability to conclude futility
                return(pmvnorm(lower = c(lk[1:(k-1)],uk[k],-Inf),
                               upper = c(uk[1:(k-1)],Inf,ck),
                               mean=c(thetheta[1:k],delta*sqrt(Info.d[k])),
                               sigma= sigmaZk2,
                               abseps = abseps) +
                       pmvnorm(lower = c(lk[1:(k-1)],-Inf,-Inf),
                               upper = c(uk[1:(k-1)],x,ck),
                               mean=c(thetheta[1:k],delta*sqrt(Info.d[k])),
                               sigma= sigmaZk2,
                               abseps = abseps) - betaSpentInc[k])
              
            }

            lk[k] <- uniroot(find.lkk,lower=uk[k]-10,upper=uk[k])$root
        }
        if(type.k %in% c("interim","decision")){
            ck <- Method1(uk=uk[1:k],
                          lk=lk[1:k],
                          Info.i=Info.i[1:k],
                          Info.d=Info.d[k],
                          Info.max=Info.max,
                          cMin=cMin,
                          ImaxAnticipated=ImaxAnticipated,
                          rho_alpha=rho_alpha,
                          alpha=alpha,
                          bindingFutility = binding)
        }else if(type.k=="final"){
            alphaSpent[k] <- alpha
            alphaSpentInc[k] <- alphaSpent[k] - alphaSpent[(k-1)]   
            betaSpent[k] <- beta
            betaSpentInc[k] <- betaSpent[k] - betaSpent[(k-1)]   
           
            uk[k] <- uniroot(function(x){pmvnorm(lower = c(TheLowerValues,x),
                                                 upper = c(uk[1:(k-1)],Inf),
                                                 mean=rep(0,k),
                                                 sigma= sigmaZk[1:k,1:k],
                                                 abseps = abseps) - betaSpentInc[k]},
                             lower = lk[utils::tail(intersect(which(!is.infinite(lk)),1:(k-1)),1)],  ## last boundary among the k-1 already computed that is not infinite  
                             upper = uk[utils::tail(intersect(which(!is.infinite(uk)),1:(k-1)),1)], 
                             tol = abseps)$root

            ## lk[k] <- uniroot(function(x){pmvnorm(lower = c(lk[1:(k-1)],-Inf),
            ##                                      upper = c(uk[1:(k-1)],x),
            ##                                      mean=thetheta[1:k],
            ##                                      sigma= sigmaZk[1:k,1:k],
            ##                                      abseps = abseps) - betaSpentInc[k]},
            ##                  lower = lk[utils::tail(intersect(which(!is.infinite(lk)),1:(k-1)),1)],  ## last boundary among the k-1 already computed that is not infinite  
            ##                  upper = uk[utils::tail(intersect(which(!is.infinite(uk)),1:(k-1)),1)],  
            ##                  tol = abseps)$root
            
            lk[k] <- uk[k]

            ck <- NA
                
            
        }
    }

    return(list(uk=uk,
                lk=lk,
                ck=ck,
                alphaSpent=alphaSpent,
                betaSpent=betaSpent))
        
}
