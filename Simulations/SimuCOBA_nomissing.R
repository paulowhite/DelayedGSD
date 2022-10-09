rm(list=ls())

name <- "Power_2_analyses_nomissing" # To save the results
path.res <- "M:\\Research\\DelayedGSD\\Github\\DelayedGSD\\Simulations\\COBA\\"
method <- 1:3 # methods used to compute the boundaries
#---
myseed <- 140786598
#--- to plan the trial ----
kMax <- 2  #max number of analyses (including final)
alpha <- 0.025  #type I error (one sided)
beta <- 0.2  #type II error
informationRates <- c(0.56,1)  #planned  information rates
rho_alpha <- 2  # rho parameter for alpha error spending function
rho_beta <- 2  # rho parameter for beta error spending function
## deltaPower <- 0.75 # just to try another value when Id > Imax
Id <- 0.67  #(expected) information rate at each decision analysis
binding <- TRUE
#
#---- to generate data -----------
#
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
PropForInterim <- 0.5 # Decide to have interim analysiz when PropForInterim % of all subjects have had the chance to have one follow-up measuement recorded in the data to be available for analysis.
theDelta.t <- 1.50001 # time lag to process the data and make them ready to analyze after collecting them (unit is time between two follow-up visits)
TimeFactor <- 14 ## number of days between two visits
#
#--- actually for both planing the trial  and generating data-----
#
#
deltaPower <- delta[3] # effect (NOT Z-scale/unit, but outcome scale/unit!) that the study is powered for: should we choose ourselves or compute from other numbers above ???
sdPower <- allsd[3]
n <- ceiling(2*2*((sdPower/deltaPower)^2)*(qnorm(1-beta)-qnorm(alpha))^2) #104 with Corine's data # should we choose ourselves or compute from the above numbers ???
# inflate SS as required for interim

#adjust for expected withdrawal
n <- n/(1-(Miss11+Miss21))

set.seed(140786598)
nsim=10000
allseeds <- sample.int(n = 1000000000, size = nsim, replace=FALSE) #x=1:(.Machine$integer.max) seems to be maximal possible

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

plannedB <- vector(mode = "list", length = 3)
for(iMeth in method){ ## iMeth <- 1
  plannedB[[iMeth]] <- CalcBoundaries(kMax=kMax,  
                                      alpha=alpha, 
                                      beta=beta,  
                                      InfoR.i=informationRates,  
                                      InfoR.d=c(Id,1),  
                                      rho_alpha=rho_alpha,  
                                      rho_beta=rho_beta,  
                                      method=iMeth,  
                                      cNotBelowFixedc=FALSE,
                                      bindingFutility=binding,
                                      delta=deltaPower)
  ## summary(plannedB[[1]])
  ## coef(plannedB[[iMeth]], type = "information")
}
inflationFactor <- unlist(lapply(plannedB,function(iP){iP$planned$InflationFactor}))
nGSD <- ceiling(n*inflationFactor)
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
                 DigitsTime=0
  )
  d <- res$d
  ## head(d,n=20)
  # }}}
  # {{{ reformat data like those of Corine
  ## Make data long format
  ## dd <- FormatAsCase(d)
  ## head(dd)
  ## summary(dd)
  # }}}
  
  # {{{ make data available at interim
  # Here we stop inclusion data collection for the interim analysis as soon as
  # half of the participants have completed (or had the opportunity to complete) the follow-up 
  thets <- d$t3[ceiling(nGSD*PropForInterim)]
  #thet <- d$t3[ceiling(n*PropForInterim)]
  
  ## ddi <- FormatAsCase(di) # needed ????
  ## head(d[d$id==52,])
  # }}}
  
  nX1.interim <- vector()
  nX2.interim <- vector()
  nX3.interim <- vector()
  currentGSD <- vector(mode = "list", length = 3)
  out.interim <- vector(mode = "list", length = 3)
  for(iMeth in method){ ## iMeth <- 1
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
    
    iConfint.interim <- confint(currentGSD[[iMeth]])
    iInfo.interim <- coef(currentGSD[[iMeth]], type = "information")
    iBoundary.interim <- coef(currentGSD[[iMeth]], type = "boundary")
    iDecision.interim <- coef(currentGSD[[iMeth]], type = "decision")
    
    out.interim[[iMeth]] <-  data.frame(statistic = iConfint.interim[1,"statistic"],
                                        estimate_ML = iConfint.interim[1,"estimate"],
                                        se_ML = iConfint.interim[1,"se"],
                                        info = iInfo.interim[1,"Interim"],
                                        infoPC = iInfo.interim[1,"Interim.pc"],
                                        info.pred = iInfo.interim[1,"Decision"],
                                        infoPC.pred = iInfo.interim[1,"Decision.pc"],
                                        uk = iBoundary.interim[1,"Ebound"],
                                        lk = iBoundary.interim[1,"Fbound"],
                                        decision = iDecision.interim["decision","stage 1"],
                                        reason = iDecision.interim["comment","stage 1"])
  }
  ## currentGSD[[1]]
  ## plot(currentGSD[[1]])
  
  
  out.decision <- vector(mode = "list", length = 3)
  for(iMeth in method){ ## iMeth <- 1
    
    ## ** decision
    dDecision <- d[which(d$t1 <= thets[iMeth] + theDelta.t*TimeFactor),]
    lmmD <- analyzeData(dDecision, ddf = "nlme", getinfo = TRUE, trace = TRUE)
    
    
    if(out.interim[[iMeth]]$decision == "stop"){
      currentGSD[[iMeth]] <- update(currentGSD[[iMeth]], delta = lmmD, trace = FALSE)
      ## plot(currentGSD[[iMeth]])
      
      iConfint.decision <- confint(currentGSD[[iMeth]], method = c("ML","MUE"))
      iInfo.decision <- coef(currentGSD[[iMeth]], type = "information")
      iBoundary.decision <- coef(currentGSD[[iMeth]], type = "boundary")
      iDecision.decision  <- coef(currentGSD[[iMeth]], type = "decision")
      
      out.decision[[iMeth]] <- data.frame(statistic = iConfint.decision[1,"statistic"],
                                          p.value_ML = iConfint.decision[iConfint.decision$method == "ML","p.value"],
                                          lower_ML = iConfint.decision[iConfint.decision$method == "ML","lower"],
                                          upper_ML = iConfint.decision[iConfint.decision$method == "ML","upper"],
                                          estimate_ML = iConfint.decision[iConfint.decision$method == "ML","estimate"],
                                          p.value_MUE = iConfint.decision[iConfint.decision$method == "MUE","p.value"],
                                          lower_MUE = iConfint.decision[iConfint.decision$method == "MUE","lower"],
                                          upper_MUE = iConfint.decision[iConfint.decision$method == "MUE","upper"],
                                          estimate_MUE = iConfint.decision[iConfint.decision$method == "MUE","estimate"],
                                          info = iInfo.decision[1,"Decision"],
                                          infoPC = iInfo.decision[1,"Decision.pc"],
                                          ck = iBoundary.decision[1,"Cbound"],
                                          decision = unname(iDecision.decision["decision","stage 1 decision"]),
                                          reason = unname(iDecision.decision["comment","stage 1 decision"]))
            
    }else{
      ## update information
      currentGSD[[iMeth]] <- update(currentGSD[[iMeth]], delta = lmmD, k = 1, type.k = "decision", trace = FALSE)
      
      iInfo.decision <- coef(currentGSD[[iMeth]], type = "information")
      iBoundary.decision <- coef(currentGSD[[iMeth]], type = "boundary")
      
      out.decision[[iMeth]] <- data.frame(statistic = NA,
                                          p.value_ML = NA,
                                          lower_ML = NA,
                                          upper_ML = NA,
                                          estimate_ML = NA,
                                          p.value_MUE = NA,
                                          lower_MUE = NA,
                                          upper_MUE = NA,
                                          estimate_MUE = NA,
                                          info = iInfo.decision[1,"Decision"],
                                          infoPC = iInfo.decision[1,"Decision.pc"],
                                          ck = iBoundary.decision[1,"Cbound"],
                                          decision = NA,
                                          reason = NA)
    }
  }
  # }}}
  # {{{ Analyze data at decision
  
  ## ** finale
  out.final <- vector(mode = "list", length = 3)
  for(iMeth in method){ ## iMeth <- 1
      dFinal <- d[1:nGSD[iMeth],]
      lmmF <- analyzeData(dFinal, ddf = "nlme", getinfo = TRUE, trace = TRUE)
    
      if(out.interim[[iMeth]]$decision == "stop"){
          out.final[[iMeth]] <- data.frame(statistic = NA,
                                           p.value_ML = NA,
                                           lower_ML = NA,
                                           upper_ML = NA,
                                           estimate_ML = NA,
                                           p.value_MUE = NA,
                                           lower_MUE = NA,
                                           upper_MUE = NA,
                                           estimate_MUE = NA,
                                           info = NA,
                                           infoPC = NA,
                                           ck = NA,
                                           decision = NA,
                                           reason = NA)
      
    }else{
      currentGSD[[iMeth]] <- update(currentGSD[[iMeth]], delta = lmmF, trace = FALSE)
      
      ## plot(test)
      ## summary(test)
      iConfint.final <- confint(currentGSD[[iMeth]], method = c("ML","MUE"))
      iInfo.final <- coef(currentGSD[[iMeth]], type = "information")
      iBoundary.final <- coef(currentGSD[[iMeth]], type = "boundary")
      iDecision.final  <- coef(currentGSD[[iMeth]], type = "decision")

      out.final[[iMeth]] <- data.frame(statistic = iConfint.final[1,"statistic"],
                                       p.value_ML = iConfint.final[iConfint.final$method == "ML","p.value"],
                                       lower_ML = iConfint.final[iConfint.final$method == "ML","lower"],
                                       upper_ML = iConfint.final[iConfint.final$method == "ML","upper"],
                                       estimate_ML = iConfint.final[iConfint.final$method == "ML","estimate"],
                                       p.value_MUE = iConfint.final[iConfint.final$method == "MUE","p.value"],
                                       lower_MUE = iConfint.final[iConfint.final$method == "MUE","lower"],
                                       upper_MUE = iConfint.final[iConfint.final$method == "MUE","upper"],
                                       estimate_MUE = iConfint.final[iConfint.final$method == "MUE","estimate"],
                                       info = iInfo.final[2,"Interim"],  #COBA: shouldn't this be taken from row 2?
                                       infoPC = iInfo.final[2,"Interim.pc"], #COBA: shouldn't this be taken from row 2?
                                       ck = iBoundary.final[2,"Cbound"], #COBA: shouldn't this be taken from row 2?
                                       decision = unname(iDecision.final["decision","stage 2"]),
                                       reason = unname(iDecision.final["comment","stage 2"])
      )
    }
  }
  # }}}
  
  stopComp <- Sys.time()
  # {{{ Save results
  
  outMerge <- do.call(rbind,lapply(method, function(iMeth){
    iNames <- unique(c(names(out.interim[[iMeth]]),names(out.decision[[iMeth]]),names(out.final[[iMeth]])))
    iMerge <- data.frame(matrix(NA, ncol = length(iNames)+3, nrow = 3, dimnames = list(NULL, c("method", "stage", "type", iNames))))
    iMerge[1,c("method","stage","type",names(out.interim[[iMeth]]))] <- data.frame(method = iMeth, stage = 1, type = "interim", out.interim[[iMeth]]) 
    iMerge[2,c("method","stage","type",names(out.decision[[iMeth]]))] <- data.frame(method = iMeth, stage = 1, type = "decision", out.decision[[iMeth]]) 
    iMerge[3,c("method","stage","type",names(out.final[[iMeth]]))] <- data.frame(method = iMeth, stage = 2, type = "final", out.final[[iMeth]])
    return(iMerge)
  }))
  
  ## outMerge[outMerge$method==3,]
  
  out <- cbind(
    ## results
    outMerge,
    ## simulation details
    time.interim = rep(thets,each=3),
    seed=myseedi,             
    nX1.interim = rep(nX1.interim,each=3),
    nX2.interim = rep(nX2.interim,each=3),
    nX3.interim = rep(nX3.interim,each=3),
    ## computation time
    computation.time=as.double(round(difftime(stopComp,startComp,units="secs"),3))
  )
  ## names(out) <- myColNames
  RES <- rbind(RES,out)
  #save(RES,file=paste0(path.res,name,"(tempo)-",nsim,".rda"))
  # }}}
}
rownames(RES) <- NULL
save(RES,file=paste0(path.res,name,"-",nsim,".rda"))

## * Summary results
load(paste0(path.res,name,"-",nsim,".rda"))
sessionInfo()
summary(RES)

res <- RES
res$final.efficacy <- res$statistic >= res$uk | res$statistic >= res$ck
res$final.efficacy[res$type%in%"interim"] <- NA

res$final.futility <- res$statistic <= res$lk | res$statistic < res$ck
res$final.futility[res$type%in%"interim"] <- NA

nsim <- length(unique(res$seed))
true_eff <- 0.6
result <- NULL
for(m in unique(res$method)){
  res_temp <- res[res$method%in%m,]
  result <- rbind(result,c("method"=m,
                           "power_bnds"=sum(res_temp$final.efficacy,na.rm=T)/nsim,
                           "power_pval"=sum(res_temp$p.value_MUE<0.025,na.rm=T)/nsim,
                           "discrep_pval_bnds"=(sum(res_temp$p.value_MUE<0.025 & !res_temp$final.efficacy,na.rm=T)+sum(!res_temp$p.value_MUE<0.05 & res_temp$final.efficacy,na.rm=T))/nsim,
                           "bias_MLE"=sum(res_temp$estimate_ML[res_temp$type%in%c("decision","final")],na.rm=T)/nsim-true_eff,
                           "bias_MUE"=sum(res_temp$estimate_MUE[res_temp$type%in%c("decision","final")]>true_eff,na.rm=T)/nsim,
                           "CI_coverage"=sum(true_eff>=res_temp$lower_MUE[res_temp$type%in%c("decision","final")] & true_eff<=res_temp$upper_MUE[res_temp$type%in%c("decision","final")],na.rm=T)/nsim
  ))
  #hist(res_temp$p.value_MUE)
}

discrep <- res[res$p.value_MUE>0.025 & !is.na(res$p.value_MUE),]
discrep <- discrep[!is.na(discrep$final.efficacy),]
dim(discrep)

discrep <- res[res$p.value_MUE<0.025 & !is.na(res$p.value_MUE),]
discrep <- discrep[is.na(discrep$final.efficacy),]
dim(discrep)



#mean(res[res$method%in%1 & res$type%in%"interim","infoPC"])
#mean(res[res$method%in%1 & res$type%in%"decision","infoPC"])
#y <- cbind(res[res$method%in%1 & res$type%in%"interim",c("info","seed")],res[res$method%in%1 & res$type%in%"decision",c("info","seed")])
#y[which(y[,1]>y[,3]),]

