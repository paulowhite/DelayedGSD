### debug-simuMain.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  9 2022 (13:51) 
## Version: 
## Last-Updated: May 12 2022 (11:13) 
##           By: Brice Ozenne
##     Update #: 8
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Method 1
library(DelayedGSD)

i.plannedB <- CalcBoundaries(kMax=2,
                             alpha=0.025,
                             beta=0.2,
                             InfoR.i=c(0.5,1.0),
                             InfoR.d=c(0.55,1),
                             rho_alpha=2,
                             rho_beta=2,
                             method=1,
                             cNotBelowFixedc=FALSE,
                             bindingFutility=TRUE,
                             delta=0.8)

i.res <- GenData(n=283,
                 N.fw=2,
                 rand.block=c(1,1,0,0),
                 allsd=c(2.5, 2.1, 2.4),
                 mean0=c(10, 0, 0),
                 delta=c(0, 0.6, 0.8),
                 ar=3.44,
                 cor.01.1=-0.15,
                 cor.ij.1=0.68,
                 cor.0j.1=-0.27,
                 seed=19519,
                 MissProb=matrix(c(0.04807692, 0.05769231, 0.00961538, 0.88461538), nrow = 2, ncol = 2,
                                 dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing"))),
                 DigitsOutcome=2,
                 TimeFactor=14,
                 DigitsTime=0
                 )

i.d <- i.res$d

## interim 
i.thet <- i.d$t3[142]
i.di <- SelectData(i.d,t=i.thet)
i.lmmI <- analyzeData(i.di, ddf = "nlme", data.decision = sum(i.d$t1 <= i.thet + 1.50001*14), getinfo = TRUE, trace = TRUE)
i.currentGSD <- update(i.plannedB, delta = i.lmmI, trace = FALSE)
i.currentGSD

## final
i.dFinal <- i.d
i.lmmF <- analyzeData(i.dFinal, ddf = "nlme", getinfo = TRUE, trace = TRUE)
i.currentGSD <- update(i.currentGSD, delta = i.lmmF, trace = FALSE)



i.currentGSD <- update.delayedGSD(i.currentGSD, delta = i.lmmF, trace = FALSE)

i.fff <- function(x){pmvnorm(lower = c(TheLowerValues,x),
                                     upper = c(uk[1:(k-1)],Inf),
                                     mean=rep(0,k),
                                     sigma= sigmaZk[1:k,1:k],
                                     abseps = abseps) - betaSpentInc[k]}

i.fff <- function(x){pmvnorm(lower = c(0.8144777,x),
                             upper = c(2.345527,Inf),
                             mean=rep(0,2),
                             sigma= cbind(c(1, 0.77139186), c(0.77139186, 1)),
                             abseps = 1e-06) - 0.1239992}

uniroot(i.fff,
        lower = 0.8144777,  ## last boundary among the k-1 already computed that is not infinite  
        upper = 2.345527, 
        tol = 1e-6)$root

uniroot(i.fff,
        lower = 0,  ## last boundary among the k-1 already computed that is not infinite  
        upper = 2.345527, 
        tol = 1e-6)$root

seqX <- seq(0,5,0.01)
plot(seqX,sapply(seqX, i.fff)) 

## * Method 3
plannedB <- CalcBoundaries(kMax=2,
                           alpha=0.025,
                           beta=0.2,
                           InfoR.i=c(0.5,1.0),
                           InfoR.d=c(0.55,1),
                           rho_alpha=2,
                           rho_beta=2,
                           method=3,
                           cNotBelowFixedc=FALSE,
                           bindingFutility=TRUE,
                           delta=0.8)

res <- GenData(n=283,
               N.fw=2,
               rand.block=c(1,1,0,0),
               allsd=c(2.5, 2.1, 2.4),
               mean0=c(10, 0, 0),
               delta=c(0, 0.6, 0.8),
               ar=3.44,
               cor.01.1=-0.15,
               cor.ij.1=0.68,
               cor.0j.1=-0.27,
               seed=19764,
               MissProb=matrix(c(0.04807692, 0.05769231, 0.00961538, 0.88461538), nrow = 2, ncol = 2,
                               dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing"))),
               DigitsOutcome=2,
               TimeFactor=14,
               DigitsTime=0
               )

d <- res$d
thet <- d$t3[142]
di <- SelectData(d,t=thet)
lmmI <- analyzeData(di, ddf = "nlme", data.decision = sum(d$t1 <= thet + 1.50001*14), getinfo = TRUE, trace = TRUE)

currentGSD <- update(plannedB, delta = lmmI, trace = FALSE)
dDecision <- d[which(d$t1 <= thet + 1.50001*14),]
lmmD <- analyzeData(dDecision, ddf = "nlme", getinfo = TRUE, trace = TRUE)

lmmD$information-lmmI$information ## decreasing information

currentGSD <- update(currentGSD, delta = lmmD, k = 1, type.k = "decision", trace = FALSE)

c_new <- uniroot(function(x){pmvnorm(lower = c(2.430755,x),
                                     upper = c(Inf,Inf),
                                     mean=rep(0,2),
                                     sigma= cbind(c(1, 1.02123651), c(1.02123651, 1)),
                                     abseps = 1e-6) - 0.007487805},
                 lower = 0.5600974,
                 upper = 2.430755,
                 tol = 1e-06)$root
##----------------------------------------------------------------------
### debug-simuMain.R ends here
