rm(list=ls())
                                   # {{{ parameters
## * parameters
NMC <- 250 # number of sequential simulations to run in parallel. Eg. with 250, then we can run 40 scripts in paralell to get N=10,000 runs in total.
kMax <- 2  #max number of analyses (including final)
alpha <- 0.025  #type I error (one sided)
beta <- 0.2  #type II error
informationRates <- c(0.5,1)  #planned  information rates
rho_alpha <- 2  # rho parameter for alpha error spending function
rho_beta <- 2  # rho parameter for beta error spending function
## deltaPower <- 0.75 # just to try another value when Id > Imax
Id <- 0.55  #(expected) information rate at each decision analysis
binding <- TRUE

block <- c(1,1,0,0) 
allsd <- c(2.5,2.1,2.4) # sd, first from baseline measurement, then the two changes from baseline
mean0 <- c(10,0,0) # mean placebo group (again, first is absolute value, then change from baseline)
delta <- c(0,0.6,0.8) # treatment effect
ar <- (0.86*2)*2 # orginial accrual rate from data from Corine is 0.86 per week, hence we multiply by 2 for by 14 days. As to low, we further multiply by 2
cor011 <- -0.15 # ~ from data from Corine
corij1 <- 0.68  # ~ from data from Corine
cor0j1 <- -0.27  # ~ from data from Corine
Miss11 <- 5/104 # miss both V1 and V2
Miss12 <- 1/104 # miss V1 and but not V2
Miss21 <- 6/104 # do not miss V1 and but miss V2
Miss22 <- 92/104 # miss none
PropForInterim <- 0.5 # Decide to have interim analysiz when PropForInterim % of all subjects have had the chance to have one follow-up measuement recorded in the data to be available for analysis.
theDelta.t <- 1.50001 # time lag to process the data and make them ready to analyze after collecting them (unit is time between two follow-up visits)
TimeFactor <- 14 ## number of days between two visits
                                        #
                                        #--- actually for both planing the trial  and generating data-----
                                        #
                                        #
deltaPower <- abs(delta[3]) # effect (NOT Z-scale/unit, but outcome scale/unit!) that the study is powered for: should we choose ourselves or compute from other numbers above ???
n <- ceiling(2*2*((allsd[3]/deltaPower)^2)*(qnorm(1-beta)-qnorm(alpha))^2) #104 with Corine's data # should we choose ourselves or compute from the above numbers ???
                                        # inflate SS as required for interim

## * libraries and functions
library(DelayedGSD)  

## * Compute inflation factor and sample size
plannedB <- CalcBoundaries(kMax=kMax,  
                           alpha=alpha, 
                           beta=beta,  
                           InfoR.i=informationRates,  
                           InfoR.d=c(Id,1),  
                           rho_alpha=rho_alpha,  
                           rho_beta=rho_beta,  
                           method=1,  
                           cNotBelowFixedc=FALSE,
                           bindingFutility=binding,
                           delta=tail(delta,1))
inflationFactor <- plannedB$planned$InflationFactor
nGSD <- ceiling(n*inflationFactor)

## * Loop
MyMissProb <- matrix(c(Miss11,Miss12,Miss21,Miss22),ncol=2,nrow=2,byrow=TRUE) # to additionnally remove 1 more because some FASFL=N
colnames(MyMissProb) <- c("V2 missing","V2 not missing")
rownames(MyMissProb) <- c("V1 missing","V1 not missing")

## ** simulate
res <- GenData(n=n, 
               N.fw=2,
               rand.block=block,
               allsd=allsd,
               mean0=mean0,
               delta=delta,
               ar=ar,
               cor.01.1=cor011,
               cor.ij.1=corij1,
               cor.0j.1=cor0j1,
               seed=8253,
               MissProb=MyMissProb,
               DigitsOutcome=2,
               TimeFactor=TimeFactor,
               DigitsTime=0
               )
d <- res$d
thet <- d$t3[ceiling(n*PropForInterim)]
di <- SelectData(d,t=thet)

## ** interim
lmmI <- analyzeData(di, ddf = "nlme", data.decision = sum(d$t1 <= thet + theDelta.t*TimeFactor), getinfo = TRUE, trace = TRUE)

currentGSD <- update(plannedB, delta = lmmI, trace = FALSE)


## ** decision
dDecision <- d[which(d$t1 <= thet + theDelta.t*TimeFactor),]
lmmD <- analyzeData(dDecision, ddf = "nlme", getinfo = TRUE, trace = TRUE)
    
## ** finale
dFinal <- d
lmmF <- analyzeData(dFinal, ddf = "nlme", getinfo = TRUE, trace = TRUE)

out.final <- vector(mode = "list", length = 3)
currentGSD <- update(currentGSD, delta = lmmF, trace = FALSE)
## Fejl i uniroot(function(x) { : 
##   f() values at end points not of opposite sign

futilityBound <- 1.04524
efficacyBound <- 2.26781
correlation <- 0.83868885
typeIIerror <- 0.1066379

uniroot(function(x){pmvnorm(lower = c(futilityBound,x),
                            upper = c(efficacyBound,Inf),
                            mean=rep(0,2),
                            sigma= cbind(c(1, correlation), 
                                         c(correlation, 1)
                                         )) - typeIIerror},
        lower = futilityBound,
        upper = efficacyBound,
        tol = 1e-6)$root

uniroot(function(x){pmvnorm(lower = c(futilityBound,x),
                            upper = c(efficacyBound,Inf),
                            mean=rep(0,2),
                            sigma= cbind(c(1, correlation), 
                                         c(correlation, 1)
                                         )) - typeIIerror},
        lower = -10,
        upper = 10,
        tol = 1e-6)$root


