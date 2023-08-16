rm(list=ls())

#load R functions
#library(devtools)
#install_github("PauloWhite/DelayedGSD")
#library(DelayedGSD)
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("R")

source("Simulations/brice-biostat-cluster/FCT.R")

#parameters
method <- 1:3 # methods used to compute the boundaries

#--- to plan the trial ----
kMax <- 2  #max number of analyses (including final)
alpha <- 0.025  #type I error (one sided)
beta <- 0.2  #type II error
informationRates <- 0.54  #planned  information rates
rho_alpha <- 2  # rho parameter for alpha error spending function
rho_beta <- 2  # rho parameter for beta error spending function
deltaPower <- 1.5 # effect (NOT Z-scale/unit, but outcome scale/unit!) that the study is powered for
sdPower <- 2.5 #sd of main outcome, used for power calculation
Id <- c(0.62,1)  #(expected) information rate at each decision analysis
binding <- FALSE
cNotBelowFixedc <- TRUE
#
#---- to generate data -----------
#
block <- c(1,1,0,0) 
allsd <- c(2.5,2.1,2.4) # sd, first from baseline measurement, then the two changes from baseline
mean0 <- c(10,0,0) # mean placebo group (again, first is absolute value, then change from baseline)
delta <- c(0,0.5,0.5) # treatment effect
delta <- c(0,0.2,0.2) # treatment effect
ar <- (0.86*2)*2 # orginial accrual rate from data from Corine is 0.86 per week, hence we multiply by 2 for by 14 days. As too low, we further multiply by 2
cor011 <- -0.15 # ~ from data from Corine
corij1 <- 0.68  # ~ from data from Corine
cor0j1 <- -0.27  # ~ from data from Corine
Miss11 <- 5/104 # miss both V1 and V2
Miss12 <- 1/104 # miss V1 and but not V2
Miss21 <- 6/104 # do not miss V1 and but miss V2
Miss22 <- 92/104 # miss none
MyMissProb <- matrix(c(Miss11,Miss12,Miss21,Miss22),ncol=2,nrow=2,byrow=TRUE, # to additionnally remove 1 more because some FASFL=N
                     dimnames = list(c("V1 missing","V1 not missing"), c("V2 missing","V2 not missing")))
PropForInterim <- 0.5 # Decide to have interim analysis when PropForInterim % of all subjects have had the chance to have one follow-up measurement recorded in the data to be available for analysis.
theDelta.t <- 1.50001 # time lag to process the data and make them ready to analyze after collecting them (unit is time between two follow-up visits)
TimeFactor <- 14 ## number of days between two visits
#
#sdPower <- sdPower*sqrt(1-cor0j1^2) #expected sd adjusted for correlation with baseline
n <- ceiling(2*2*((sdPower/deltaPower)^2)*(qnorm(1-beta)-qnorm(alpha))^2) #104 with Corine's data # should we choose ourselves or compute from the above numbers ???
# inflate SS as required for interim

#adjust for expected withdrawal
#n <- n/(1-(Miss11+Miss21))
n <- n/(1-0.18) #54/arm based on planned withdrawal rate

#if we do an interim after 50% of the patients has (had the chance to) complete the Day 28 visit, this means
ar*2
#more patients will be recruited while waiting for the last patient to complete. Half of these patients will contribute with at least one visit
#a rough estimate of the information fraction at interim could be
(n/2 + ar)/n
#so around 0.56 

#the expected information at the decision analysis is around the following information fraction
(theDelta.t*ar+ar*2+n/2)/n  #half of the patients recruited, plus those recruited during 28 days plus those recruited during 3 weeks data processing
#0.62

#planned boundaries
plannedB <- vector(mode = "list", length = 3)
for(iMeth in method){ ## iMeth <- 1
  plannedB[[iMeth]] <- CalcBoundaries(kMax=kMax,  
                                      alpha=alpha, 
                                      beta=beta,  
                                      InfoR.i=informationRates,  
                                      InfoR.d=Id,  
                                      rho_alpha=rho_alpha,  
                                      rho_beta=rho_beta,  
                                      method=iMeth,  
                                      cNotBelowFixedc=cNotBelowFixedc,
                                      bindingFutility=binding,
                                      delta=deltaPower)
  ## summary(plannedB[[1]])
  ## coef(plannedB[[iMeth]], type = "information")
}
inflationFactor <- unlist(lapply(plannedB,function(iP){iP$planned$InflationFactor}))
nGSD <- ceiling(n*inflationFactor)

#correcting as in example t-test power calculation used, plus equal number per arm so cannot have uneven total sample size
nGSD <- c(116,116,116)

png("PlannedBndsExamplePaper.png")
plot(plannedB[[1]])
dev.off()

#generate data:
myseedi <- 56922#34513

res <- GenData(n=max(nGSD), 
               N.fw=2,
               rand.block=block,
               allsd=allsd,
               mean0=mean0,
               delta=delta,
               ar=ar,
               cor.01.1=cor011,
               cor.ij.1=corij1,
               cor.0j.1=cor0j1,
               seed=myseedi,
               MissProb=MyMissProb,
               DigitsOutcome=2,
               TimeFactor=TimeFactor,
               DigitsTime=0)

d <- res$d

# Here we stop inclusion data collection for the interim analysis as soon as
# half of the participants have completed (or had the opportunity to complete) the follow-up 
thets <- d$t3[ceiling(nGSD*PropForInterim)]

nX1.interim <- vector()
nX2.interim <- vector()
nX3.interim <- vector()
currentGSD <- vector(mode = "list", length = 3)
out.interim <- vector(mode = "list", length = 3)
for(iMeth in method){ ## iMeth <- 1
  #iMeth <- 1
  # {{{ make data available at interim
  di <- SelectData(d,t=thets[iMeth])
  
  nX1.interim[iMeth] <- sum(!is.na(di$X1))
  nX2.interim[iMeth] = sum(!is.na(di$X2))
  nX3.interim[iMeth] = sum(!is.na(di$X3))
  
  ## {{{ analyze data at at interim
  ## ** interim
  lmmI <- analyzeData(di, ddf = "nlme", data.decision = sum(d$t1 <= thets[iMeth] + theDelta.t*TimeFactor), getinfo = TRUE, trace = TRUE)
  ## lmmI <- analyzeData(di, ddf = "nlme", getinfo = TRUE, trace = TRUE)
  
  currentGSD[[iMeth]] <- update(plannedB[[iMeth]], delta = lmmI, trace = FALSE)
  
  out.interim[[iMeth]] <- exportGSD(currentGSD[[iMeth]],
                                    export.statistic = TRUE,
                                    export.ML = TRUE,
                                    export.MUE = FALSE,
                                    export.info = TRUE,
                                    export.predinfo = TRUE,
                                    export.boundary = TRUE,
                                    export.decision = TRUE)
}

png("InterimBndsExamplePaper.png")
plot(currentGSD[[1]])
dev.off()

out.decision <- vector(mode = "list", length = 3)
for(iMeth in method){ ## iMeth <- 1
  
  ## ** decision
  dDecision <- d[which(d$t1 <= (thets[iMeth] + theDelta.t*TimeFactor)),]
  lmmD <- analyzeData(dDecision, ddf = "nlme", getinfo = TRUE, trace = TRUE)
  
  if(out.interim[[iMeth]]$decision == "stop"){
    currentGSD[[iMeth]] <- update(currentGSD[[iMeth]], delta = lmmD, trace = FALSE)
    ## plot(currentGSD[[iMeth]])
    
    out.decision[[iMeth]] <- exportGSD(currentGSD[[iMeth]],
                                       export.statistic = TRUE,
                                       export.ML = TRUE,
                                       export.MUE = TRUE,
                                       export.info = TRUE,
                                       export.predinfo = FALSE,
                                       export.boundary = TRUE,
                                       export.decision = TRUE)
    
    
  }else{
    ## update information
    currentGSD[[iMeth]] <- update(currentGSD[[iMeth]], delta = lmmD, k = 1, type.k = "decision", trace = FALSE)
    
    out.decision[[iMeth]] <- exportGSD(currentGSD[[iMeth]],
                                       export.statistic = FALSE,
                                       export.ML = FALSE,
                                       export.MUE = FALSE,
                                       export.info = TRUE,
                                       export.predinfo = FALSE,
                                       export.boundary = TRUE,
                                       export.decision = FALSE)
  }
}

FinalPvalue(2.846514,  
            2.55245,  
            1.959964,
            1.739491,   
            1.085299,  
            2.261545,
            "futility",
            kMax=2, 
            delta=0,  
            estimate=0.3476515,
            method=1,  
            bindingFutility=TRUE,
            cNotBelowFixedc=TRUE,
            continuity.correction=TRUE)

FinalPvalue(2.846514,  
            2.55245,  
            1.959964,
            1.923053,   
            1.047626,  
            2.216121,
                        "futility",
                        kMax=2, 
                        delta=0,  
                        estimate=0.3476515,
                        method=3,  
                        bindingFutility=TRUE,
                        cNotBelowFixedc=TRUE,
                        continuity.correction=TRUE)

png("DecisionBndsExamplePaper.png")
plot(currentGSD[[1]])
dev.off()
