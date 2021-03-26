#Example values that can be used
#bnds <- CalcBoundaries(kMax=3,informationRates = c(0.5,0.75,1),Id=c(0.55,0.8))
#
#Id <- bnds$Id
#Ik <- bnds$Ik[1:2]
#ck <- bnds$ck
#delta <- bnds$delta
#sided=bnds$sided
#estimate=2
#lk <- bnds$lk[1:2]
#uk <- bnds$uk[1:2]
#kMax <- 3

FinalPvalue <- function(Id,  #Information at all decision analyses up to stage where study was stopped
                        Ik,  #Information at all interim analyses up to stage where study was stopped
                        ck,   #decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis)
                        lk,  #lower bounds up to stage where study was stopped
                        uk,  #upper bounds up to stage where study was stopped
                        sided=1,  #whether the test is 1 or 2-sided
                        kMax, #maximum number of analyses
                        delta=0,  #true effect under which to calculate the probability (should always be 0 for p-value, only change for calculation of CI)
                        estimate){ #the observed treatment estimate at decision
  require(mvtnorm)
  
  if(sided!=1){
    stop("Function cannot handle two-sided tests yet")
  }
  
  ## message("the method assumes that positive effects are good")
  k <- length(Id)
  I_all <- c(rbind(Ik, Id))
  
  m <- length(I_all)
  sigmaZm <- diag(1,m)
  for(i in 1:m){
    for(j in i:m){
      sigmaZm[i,j] <- sqrt(I_all[i]/I_all[j])
      sigmaZm[j,i] <- sqrt(I_all[i]/I_all[j])
    }
  }
  
  theta <- delta*sqrt(I_all)
  
  index <- which(Id > Id[k])  #find decision information levels that are higher than final information levels.
  
  if(1 %in% index){
    pval <- 0
  } else {
    pval <- pmvnorm(lower = c(uk[1],ck[1],rep(-Inf,m-2)),  #prob to stop for eff at first analysis and conclude eff
                    upper = rep(Inf,m),
                    mean=theta,
                    sigma= sigmaZm) +
      pmvnorm(lower = c(-Inf,ck[1],rep(-Inf,m-2)),   #prob to stop for fut at first analysis and switch to eff
              upper = c(lk[1],rep(Inf,m-1)),
              mean=theta,
              sigma= sigmaZm)
  }
  
  if(k>1){
    for(i in 2:k){
      
      if(i%in%index & estimate*sqrt(Id[k])>=ck[k]){  #results for larger information levels are never more extreme than smaller information levels with results above c
        next
      } else if(i %in% index){  #results for larger information levels are always more extreme if observed result < c
       #i cannot be kMax in this case, since Id[k] would be Id[kMax] in the definition of the index
         pval <- pval +  pmvnorm(lower = c(lk[1:(i-1)],uk[i],-Inf,rep(-Inf,m-1-i)),  #prob to stop for eff at analysis i and then obtain whichever result
                                upper = c(uk[1:(i-1)],rep(Inf,m-i+1)),
                                mean=theta,
                                sigma= sigmaZm) +
                         pmvnorm(lower = c(lk[1:(i-1)],-Inf,-Inf,rep(-Inf,m-1-i)),   #prob to stop for fut at analysis i and then obtain whichever result
                                upper = c(uk[1:(i-1)],lk[i],rep(Inf,m-i)),
                                mean=theta,
                                sigma= sigmaZm)
      } else {  #if information is always increasing
        if(i<kMax){
          pval <- pval +  pmvnorm(lower = c(lk[1:(i-1)],uk[i],ck[i],rep(-Inf,m-1-i)),  #prob to stop for eff at analysis i and conclude eff
                                  upper = c(uk[1:(i-1)],rep(Inf,m-i+1)),
                                  mean=theta,
                                  sigma= sigmaZm) +
                          pmvnorm(lower = c(lk[1:(i-1)],-Inf,ck[1],rep(-Inf,m-1-i)),   #prob to stop for fut at analysis i and switch to eff
                                  upper = c(uk[1:(i-1)],lk[i],rep(Inf,m-i)),
                                  mean=theta,
                                  sigma= sigmaZm)
        } else if(i==kMax){
          pval <- pval + pmvnorm(lower = c(lk[1:(i-1)],estimate*sqrt(Ik[kMax])),  #prob to continue until final analysis and obtain a higher result than the observed
                                 upper = c(uk[1:(i-1)],Inf),
                                 mean=theta,
                                 sigma= sigmaZm)
        } else {
          stop("Length of Id > kMax")
        }
      }
    }
  }
  pval
}


######try to match with rpact (almost there up to 5th decimal)
#library(rpact)
#design <- getDesignGroupSequential(kMax=3, sided = 1, alpha = 0.025, beta = 0.2,
#                         informationRates = c(0.5,0.75,1),
#                         typeOfDesign="asKD",
#                         typeBetaSpending="bsKD",gammaA=2,gammaB=2)

#res <- getDataset(sampleSizes1 = c(40,20),sampleSizes2=c(40,20),means1=c(1.5,3),means2=c(0.5,-1),stDevs1=c(2,1.65),stDevs2=c(2,1.65))
#ests <- getAnalysisResults(design,res,normalApproximation=TRUE)

#pvalue=0.006281
#CI=0.238,1.878
#estimate=1.086

#means diffs: 1,2
#SEs: sqrt(2^2/40 + 2^2/40), sqrt(2^2/60 + 2^2/60)
#test statistics 2.236, 3.162
#critvals: 2.498,2.292
#futility bnds 0.443, 1.270

#k <- 2
##Ik <- c(0.5,0.75)
#Ik <- c(1/(2^2/40 + 2^2/40), 1/(2.007^2/60 + 2.007^2/60))
#sigmaZk <- diag(1,2)
#for(i in 1:k){
#  for(j in i:k){
#    sigmaZk[i,j] <- sqrt(Ik[i]/Ik[j])
#    sigmaZk[j,i] <- sqrt(Ik[i]/Ik[j])
#  }
#}

#need to calculate the probability of a more extreme result
#that means: stopping at first analysis for efficacy
#            having an even higher result at the second analysis without stopping at the first analysis
#library(mvtnorm)

#pmvnorm(lower = c(2.498,-Inf),
#        upper = c(Inf,Inf),
#        mean=rep(0,2),
#        sigma= sigmaZk)+
#pmvnorm(lower = c(0.443,5.457),   #3.162
#        upper = c(2.498,Inf),
#        mean=rep(0,2),
#        sigma= sigmaZk)

