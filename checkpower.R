rm(list=ls())

beta <- 0.2
alpha <- 0.025

block <- c(1,1,0,0) 
allsd <- c(2.5,2.1,2.4) # sd, first from baseline measurement, then the two changes from baseline
mean0 <- c(10,0,0) # mean placebo group (again, first is absolute value, then change from baseline)
delta <- c(0,0.3,0.6) # treatment effect
ar <- (0.86*2)*2*5 # orginial accrual rate from data from Corine is 0.86 per week, hence we multiply by 2 for by 14 days. As too low, we further multiply by 2
cor011 <- -0.15 # ~ from data from Corine
corij1 <- 0.68  # ~ from data from Corine
cor0j1 <- -0.27  # ~ from data from Corine
Miss11 <- 0#5/104 # miss both V1 and V2
Miss12 <- 0#1/104 # miss V1 and but not V2
Miss21 <- 0#6/104 # do not miss V1 and but miss V2
Miss22 <- 1#92/104 # miss none
theDelta.t <- 1.50001 # time lag to process the data and make them ready to analyze after collecting them (unit is time between two follow-up visits)
TimeFactor <- 14 ## number of days between two visits


deltaPower <- delta[3] # effect (NOT Z-scale/unit, but outcome scale/unit!) that the study is powered for: should we choose ourselves or compute from other numbers above ???
sdPower <- allsd[3]
n <- ceiling(2*2*((sdPower/deltaPower)^2)*(qnorm(1-beta)-qnorm(alpha))^2) #104 with Corine's data # should we choose ourselves or compute from the above numbers ???
n <- n/(1-(Miss11+Miss21))

set.seed(140786598)
nsim=10000
allseeds <- sample.int(n = 1000000000, size = nsim, replace=FALSE) #x=1:(.Machine$integer.max) seems to be maximal possible

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("R")

RES <- NULL

## * Loop
allj <- 1:nsim
#allj <- 572:1000
for(j in allj){ ## j <- 51 ## 5
  startComp <- Sys.time()
  myseedi <- allseeds[j]
  #myseedi <- 474821724
  # {{{ TRACE info (e.g. to check the Rout)
  print(paste0("seed ",myseedi," for ","j=",which(j==allj)," out of ",nsim))
  # }}}
  # {{{ Missing probabilities
  MyMissProb <- matrix(c(Miss11,Miss12,Miss21,Miss22),ncol=2,nrow=2,byrow=TRUE) # to additionnally remove 1 more because some FASFL=N
  colnames(MyMissProb) <- c("V2 missing","V2 not missing")
  rownames(MyMissProb) <- c("V1 missing","V1 not missing")
  # }}}
  
  # {{{ generate data
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
                 seed=myseedi,
                 MissProb=MyMissProb,
                 DigitsOutcome=2,
                 TimeFactor=TimeFactor,
                 DigitsTime=0
  )
  d <- res$d
  
  lmm <- analyzeData(d, ddf = "nlme", trace = TRUE)
  RES <- c(RES,lmm$delta$p.value<alpha)
}

mean(RES)
save(RES,file="checkpowersims_10000.Rdata")
#0.8264