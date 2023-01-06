### example-mismatch-pVal-boundary.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 13 2022 (14:49) 
## Version: 
## Last-Updated: jan  6 2023 (14:46) 
##           By: Brice Ozenne
##     Update #: 14
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Simulate data
library(DelayedGSD)

method <- 3
binding <- TRUE
seed <- 65836753

debug.plannedB <- CalcBoundaries(kMax = 2,  
                                 alpha = 0.025, 
                                 beta = 0.2,  
                                 InfoR.i = c(0.482, 1.00),  
                                 InfoR.d = c(0.52,1),  
                                 rho_alpha = 2,  
                                 rho_beta = 2,  
                                 method = method,  
                                 cNotBelowFixedc = FALSE,
                                 bindingFutility= binding,
                                 delta = 0.6)

df.sim <- GenData(n = 486, 
                  N.fw = 2,
                  rand.block = c(1,1,0,0),
                  allsd = c(2.5,2.1,2.4),
                  mean0 = c(10,0,0),
                  delta = c(0,0.3,0.6),
                  ar = (0.86*2)*2*5,
                  cor.01.1 = -0.15,
                  cor.ij.1 = 0.68,
                  cor.0j.1 = -0.27,
                  seed = seed,
                  MissProb = matrix(c(0, 0, 0, 0), 
                                    nrow = 2, 
                                    ncol = 2, 
                                    dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing")) 
                                    ) ,
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
## plot(debug.GSD)
## summary(debug.GSD, planned = FALSE, predicted = FALSE)

## stage |         Interim           | Decision       |    Spent         
##       | F-bound E-bound    Stat   |  C-bound Stat  |    alpha     beta
##     1 | 0.22425 2.52429 2.13162 C |  1.95999       | 0.005447 0.043576
##     2 |                           |  1.99777       |    0.025      0.2

## Observed and planned information: 
## stage |  Interim     (%) | Decision    (%) |       n
##     1 | 10.60272 0.46678 | 12.30914 0.5419 | 299.000
##     2 |                  | 22.71478 1.0000 | 640.564

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

    df.sim2 <- df.sim[1:100,]
    df.sim2$id <- paste0("A",df.sim2$id)
    debug.lmmF <- analyzeData(rbind(df.sim[1:nGSD[method],], df.sim2),
                              ddf = "nlme", getinfo = TRUE, trace = TRUE)
    debugF.GSD <- update(debug.GSD, delta = debug.lmmF, trace = FALSE)

    debugF.GSD$ck
    
    ## alternative calculation of the p-value (Paul's version)
    term1 <- debugF.GSD$alphaSpent[1]
    term2 <- pmvnorm(lower = c(0.22425, 2.25911),
                     upper = c(2.52429353,Inf),
                     mean=c(0,0),
                     sigma= matrix(c(1,sqrt(10.60271846/24.67093),sqrt(10.60271846/24.67093),1),2,2))
    term1 + term2 - confint(debugF.GSD)[2,"p.value"]

}    


## * By hand
Info.d <- 12.58798
Info.i <- c(10.60272, 20.74867)
ck <- 1.959964
lk <- c(0.2242486, 1.9939758) 
uk <-  c(2.524294, 1.993976)  
thetaHat <- 0.4374693
alphaSpent <- debug.GSD$alphaSpent[1]

sigmaZm <- diag(1,3,3)
sigmaZm[1,2] <- sigmaZm[2,1] <- sqrt(Info.i[1]/Info.d)
sigmaZm[1,3] <- sigmaZm[3,1] <- sqrt(Info.i[1]/Info.i[2])
sigmaZm[2,3] <- sigmaZm[3,2] <- sqrt(Info.d/Info.i[2])

corrected.ck <- uniroot(function(x){
    pmvnorm(lower = c(uk[1],x),
            upper = c(Inf,Inf),
            mean=rep(0,2),
            sigma= sigmaZm[c(1,2),c(1,2)],
            abseps = 1e-6) - alphaSpent},
    lower = 0,
    upper = 3,
    tol = 1e-6)$root

    
## ** Compute boundary
uniroot(function(x){pmvnorm(lower = c(lk[1],x),
                            upper = c(uk[1],Inf),
                            mean=rep(0,2),
                            sigma= sigmaZm[c(1,3),c(1,3)],
                            abseps = 1e-6) - (0.025-alphaSpent)},
        lower = 0,
        upper = 3,
        tol = 1e-6)$root
## [1] 1.993976
pmvnorm(lower = c(lk[1],uk[2]),
        upper = c(uk[1],Inf),
        mean=rep(0,2),
        sigma= sigmaZm[c(1,3),c(1,3)],
        abseps = 1e-6)
## [1] 0.01955299

## ** Compute p-value
resP <- FinalPvalue(Info.d = Info.d,  
                    Info.i = Info.i,
                    ck = ck, 
                    lk = lk,  
                    uk = uk,  
                    kMax = 2, 
                    delta = 0,  
                    estimate = thetaHat,
                    method = 3,
                    bindingFutility = TRUE,
                    cNotBelowFixedc = FALSE)
resP < 0.025
thetaHat*sqrt(Info.i[2]) > lk[2]

FinalPvalue(Info.d = Info.d,  
            Info.i = Info.i,
            ck = corrected.ck, 
            lk = lk,  
            uk = uk,  
            kMax = 2, 
            delta = 0,  
            estimate = thetaHat,
            method = 3,
            bindingFutility = TRUE,
            cNotBelowFixedc = FALSE)

debug.GSD$alphaSpent[1] + pContEff
pStopEffEff + pContEff


pStopEffEff <- mvtnorm::pmvnorm(lower = c(uk[1],ck), 
                                upper = c(Inf,Inf),
                                mean = c(0,0),
                                sigma= sigmaZm[1:2,1:2,drop=FALSE])
pStopEffEff.original <- mvtnorm::pmvnorm(lower = c(uk[1],corrected.ck), 
                                         upper = c(Inf,Inf),
                                         mean = c(0,0),
                                         sigma= sigmaZm[1:2,1:2,drop=FALSE])
pStopEffFu <- mvtnorm::pmvnorm(lower = c(uk[1],-Inf), 
                               upper = c(Inf,ck),
                               mean = c(0,0),
                               sigma= sigmaZm[1:2,1:2,drop=FALSE])
pStopFu <- mvtnorm::pmvnorm(lower = c(-Inf,-Inf), 
                            upper = c(lk[1],Inf),
                            mean = c(0,0),
                            sigma= sigmaZm[1:2,1:2,drop=FALSE])
pContEff <- mvtnorm::pmvnorm(lower = c(lk[1], thetaHat*sqrt(Info.i[2])), 
                             upper = c(uk[1],Inf),
                             mean = c(0,0),
                             sigma= sigmaZm[c(1,3),c(1,3),drop=FALSE])
pContFu <- mvtnorm::pmvnorm(lower = c(lk[1],-Inf), 
                            upper = c(uk[1],thetaHat*sqrt(Info.i[2])),
                            mean = c(0,0),
                            sigma= sigmaZm[c(1,3),c(1,3),drop=FALSE])

pStopEffEff + pStopEffFu + pStopFu + pContEff + pContFu


pStopEffEff - alphaSpent[1]


pStopEffEff.original + pContEff ## more extreme
pStopEffEff + pContEff ## more extreme
1 - (pStopFu + pStopEffFu + pContFu) ## less extreme  

##----------------------------------------------------------------------
### example-mismatch-pVal-boundary.R ends here
