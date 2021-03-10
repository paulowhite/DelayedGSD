### SimuMain.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar  5 2021 (10:56) 
## Version: 
## Last-Updated: Mar 10 2021 (14:12) 
##           By: Paul Blanche
##     Update #: 375
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


rm(list=ls())


# {{{ parameters
name <- "ScenarioName" # To save the results
#----
NMC <- 10 # number of sequential simulations to run in parallel. Eg. with 250, then we can run 40 scripts in paralell to get N=10,000 runs in total.
#---
myseed <- 140786598
#--- to plan the trial ----
kMax <- 2  #max number of analyses (including final)
alpha <- 0.025  #type I error (one sided)
beta <- 0.2  #type II error
informationRates <- c(0.5,1)  #planned or observed information rates
gammaA <- 2  # rho parameter for alpha error spending function
gammaB <- 2  # rho parameter for beta error spending function
## deltaPower <- 0.75 # just to try another value when Id > Imax
Id <- 0.55  #(expected) information ratio at each decision analysis
#
#---- to generate data -----------
#
block <- c(1,1,0,0) 
allsd <- c(2.5,2.1,2.4)
mean0 <- c(10,0,0) # 
delta <- c(0,-0.6,-0.8) # 
ar <- (0.86*2)*2 # orginial accrual rate from data from Corine is 0.86 per week, hence we multiply by 2 for by 14 days. As to low, we further multiply by 2
cor011 <- -0.15 # ~ from data from Corine
corij1 <- 0.68  # ~ from data from Corine
cor0j1 <- -0.27  # ~ from data from Corine
Miss11 <- 5/104 # miss both V1 and V2
Miss12 <- 1/104 # miss V1 and but not V2
Miss21 <- 6/104 # do not miss V1 and but miss V2
Miss22 <- 92/104 # miss none
PropForInterim <- 0.5 # Decide to have interim analysiz when PropForInterim % of all subjects have had the chance to have one follow-up measuement recorded in the data to be available for analysis.
theDelta.t <- 0.50001 # time lag to process the data and make them ready to analyze after collecting them (unit is time between two follow-up visits)
#
#--- actually for both planing the trial  and generating data-----
#
#
deltaPower <- abs(delta[3]) # effect (NOT Z-scale/unit, but outcome scale/unit!) that the study is powered for: should we choose ourselves or compute from other numbers above ???
n <- ceiling(2*2*((allsd[3]/deltaPower)^2)*(qnorm(1-beta)-qnorm(alpha))^2) #104 with Corine's data # should we choose ourselves or compute from the above numbers ???
# inflate SS as required for interim
R <- getDesignCharacteristics(getDesignGroupSequential(kMax=kMax, sided = 1, alpha = alpha, beta = beta,
                                                       informationRates = informationRates,
                                                       typeOfDesign="asKD",
                                                       typeBetaSpending="bsKD",gammaA=gammaA,gammaB=gammaB))$inflationFactor
n <- ceiling(n*R)
#
# --- just to check---
## n
## InfoFixed <- ((qnorm(1-beta)-qnorm(alpha))/(deltaPower))^2
## InfoFixed
## n/(4*(allsd[3])^2)
#---
# }}}

# {{{ set path to load and save and othe machine specific variables
# if I am on my own laptop
if(system("whoami",intern=TRUE)=="paul"){  
    pathToLoad <- "~/research/SeqDesignDelayed/DelayedGSD/Rfunctions/"
    pathToSave <- "~/research/SeqDesignDelayed/DelayedGSD/Simulations/output/"
    name <- paste0(name,"test")
    i <- 1
}
# if I am on the biostat server
if(system("whoami",intern=TRUE)=="tfq625"){
    pathToLoad <- "~/to-sync-server-biostat/TOCREATE"
    pathToSave <- "~/to-sync-server-biostat/TOCREATE"
    i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}
# }}}

# {{{ Source helper functions
sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}
sourceDir(pathToLoad)
# }}}

# {{{ Set seeds for parallel computing and reproducibility
## RNGkind(sample.kind = "Rounding")  # to reproduce old simulation results !!!!!!!!!!! (the default random number generator has been updated with in R version 3.6)
## RNGkind(sample.kind = "Rejection") # new default
set.seed(myseed)
allseeds <- sample(x=1:1e8,size=1e5,replace=FALSE) #x=1:(.Machine$integer.max) seems to be maximal possible
# unique(round(runif(n=1e5)*1e10))
## length(allseeds)
## head(allseeds)
## sample(x=1:1e10,size=1e5,replace=FALSE)
# }}}

# {{{ technical details to loop
RES <- NULL # initialize results to save
allj <- ((i-1)*NMC + 1):(i*NMC) # indices of all iterations (replicates) for this job, accountng for the other jobs running in parallel
# }}}

# {{{ Compute planned bundries
#
#-- First case: method 1,  cutoff at decision can be c < 1.96
#
PlannedB1 <- CalcBoundaries(kMax=kMax,  #max number of analyses (including final)
                            sided=1,  #one or two-sided
                            alpha=alpha,  #type I error
                            beta=beta,  #type II error
                            informationRates=informationRates,  #planned or observed information rates
                            gammaA=gammaA,  #rho parameter for alpha error spending function
                            gammaB=gammaB,  #rho parameter for beta error spending function
                            method=1,  #use method 1 or 2 from paper H&J
                            cNotBelowFixedc=FALSE, # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                            delta=deltaPower,  #effect that the study is powered for
                            Id=Id)    #(expected) information ratio at each decision analysis
#
#-- Second case: method 1, cutoff at decision c >= 1.96
#
PlannedB2 <- CalcBoundaries(kMax=kMax,  #max number of analyses (including final)
                            sided=1,  #one or two-sided
                            alpha=alpha,  #type I error
                            beta=beta,  #type II error
                            informationRates=informationRates,  #planned or observed information rates
                            gammaA=gammaA,  #rho parameter for alpha error spending function
                            gammaB=gammaB,  #rho parameter for beta error spending function
                            method=1,  #use method 1 or 2 from paper H&J
                            cNotBelowFixedc=TRUE, # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                            delta=deltaPower,  #effect that the study is powered for
                            Id=Id)    #(expected) information ratio at each decision analysis
#---

# }}}

for(j in allj){
    startComp <- Sys.time()
    # {{{ TRACE info (e.g. to check the Rout)
    print(paste0("j=",which(j==allj)," out of ",NMC, " (i.e. j=",j," in [",(i-1)*NMC + 1,";",i*NMC,"], as job id is i=",i,")"))
    ## myseedi <- res.mat[i,"seed"]
    myseedi <- allseeds[j]
    print(paste0("seed=",myseedi))    
    # }}}
    # {{{ Missing probabilities
    MyMissProb <- matrix(c(Miss11,Miss12,Miss21,Miss22),ncol=2,nrow=2,byrow=TRUE) # to additionnally remove 1 more because some FASFL=N
    colnames(MyMissProb) <- c("V2 missing","V2 not missing")
    rownames(MyMissProb) <- c("V1 missing","V1 not missing")
    # }}}

    # {{{ generate data
    set.seed(myseedi)
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
                   TimeFactor=14,
                   DigitsTime=0
                   )
    d <-res$d
    head(d,n=20)
    # }}}
    # {{{ reformat data like those of Corine
    ## Make data long format
    dd <- FormatAsCase(d)
    head(dd)
    ## summary(dd)
    # }}}
    # {{{ make data available at interim
    thet <- d$t2[ceiling(n*PropForInterim)] + theDelta.t
    di <- SelectData(d,t=thet)
    ddi <- FormatAsCase(di) # needed ????
    ## head(d[d$id==52,])
    # }}}
    # {{{ analyze data at at interim
    ResInterim <- AnalyzeData(di) #  Log-restricted-likelihood: -189.561
    # save (potentially) important results
    Z.interim <- ResInterim$estimate/ResInterim$se
    est.interim <- ResInterim$estimate
    se.interim <- ResInterim$se
    Info.interim <- ResInterim$Info
    PrCtInfo.interim <- ResInterim$Info/PlannedB1$Imax # Information rate at interim (denominator is targetted Imax)
    ## dittes <- di[di$id <=52,] just for checking how missing data are handles for subjects with missing at both follow-up times
    ## ResInterimtes <- AnalyzeData(dittes) # OK, log-likelihood is the same
    # }}}


    # {{{ Predict info at Decision if we stop at interim, based on current data (for method 2 only)
    ## getInformation(Res$fit,
                   ## data = Res$d.long,
                   ## name.coef = "Z1",
                   ## type = "prediction",
                   ## method.prediction = "inflation")
    # }}}
    
    # {{{ Decision at Interim
    DecisionInterimB1 <- Decision(ResInterim,  #results from AnalyzeData
                                  PlannedB1,  #results from CalcBoundaries
                                  k=1, #at which phase are we?
                                  analysis="interim", #is it an interim or decision or final analysis
                                  Ik=c(ResInterim$Info,PlannedB1$Imax),  #all I_k (information) from first interim to final analysis (observed where possible)
                                  Id=NULL,  #expected or observed information at each decision analysis
                                  PositiveIsGood=FALSE, # whether positive effect is good (i.e. positive trial)
                                  Trace=FALSE, # whether to print some messages
                                  plot=FALSE)  #should the boundaries and results be plotted?
    # save (potentially) important results
    u.interim.B1 <- DecisionInterimB1$details["u"]
    l.interim.B1 <- DecisionInterimB1$details["l"]
    dec.interim.B1 <- DecisionInterimB1$decision
    #---
    DecisionInterimB2 <- Decision(ResInterim,  #results from AnalyzeData
                                  PlannedB2,  #results from CalcBoundaries
                                  k=1, #at which phase are we?
                                  analysis="interim", #is it an interim or decision or final analysis
                                  Ik=c(ResInterim$Info,PlannedB2$Imax),  #all I_k (information) from first interim to final analysis (observed where possible)
                                  Id=NULL,  #expected or observed information at each decision analysis
                                  PositiveIsGood=FALSE, # whether positive effect is good (i.e. positive trial)
                                  Trace=FALSE, # whether to print some messages
                                  plot=FALSE)  #should the boundaries and results be plotted?
    # save (potentially) important results
    u.interim.B2 <- DecisionInterimB2$details["u"]
    l.interim.B2 <- DecisionInterimB2$details["l"]
    dec.interim.B2 <- DecisionInterimB2$decision   
    # }}}
   
    # {{{ Decision at interim is to stop ("Efficacy" or "Futility"): move to decision analysis
    if(DecisionInterimB1$decision!="Continue"){
        # Make data available at Decision Analysis
        dDecision <- d[which(d$id %in% di$id),]
        # Analyze data at at interim
        ResDecision <- AnalyzeData(dDecision)   
        # Decision at Decision Analysis
        DecisionDecisionB1 <- Decision(ResDecision,  #results from AnalyzeData
                                       PlannedB1,  #results from CalcBoundaries
                                       k=1, #at which phase are we?
                                       analysis="decision", #is it an interim or decision or final analysis
                                       Ik=c(ResInterim$Info,PlannedB1$Imax),  #all I_k (information) from first interim to final analysis (observed where possible)
                                       Id=max(ResDecision$Info,ResInterim$Info+0.01)/PlannedB1$Imax,  #expected or observed information RATIO at each decision analysis
                                       PositiveIsGood=FALSE, # whether positive effect is good (i.e. positive trial)
                                       Trace=FALSE, # whether to print some messages
                                       plot=FALSE)  #should the boundaries and results be plotted?
        DecisionDecisionB2 <- Decision(ResDecision,  #results from AnalyzeData
                                       PlannedB2,  #results from CalcBoundaries
                                       k=1, #at which phase are we?
                                       analysis="decision", #is it an interim or decision or final analysis
                                       Ik=c(ResInterim$Info,PlannedB2$Imax),  #all I_k (information) from first interim to final analysis (observed where possible)
                                       Id=max(ResDecision$Info,ResInterim$Info+0.01)/PlannedB2$Imax,  #expected or observed information RATIO at each decision analysis
                                       PositiveIsGood=FALSE, # whether positive effect is good (i.e. positive trial)
                                       Trace=FALSE, # whether to print some messages
                                       plot=FALSE)  #should the boundaries and results be plotted?
        # to save
        ConclusionTrial.B1 <- DecisionDecisionB1$decision
        c.decision.B1 <- DecisionDecisionB1$details
        ConclusionTrial.B2 <- DecisionDecisionB2$decision
        c.decision.B2 <- DecisionDecisionB2$details
        Z.decision <- ResDecision$estimate/ResDecision$se
        est.decision <- ResDecision$estimate
        se.decision <- ResDecision$se
        Info.decision <- ResDecision$Info
        PrCtInfo.decision  <-  ResDecision$Info/PlannedB1$Imax
    }else{
        c.decision.B1 <- NA
        c.decision.B2 <- NA
        Z.decision <- NA
        est.decision <- NA
        se.decision <- NA
        Info.decision <- NA
        PrCtInfo.decision  <-  NA

    }
    # }}}

    # {{{ Decision is to continue: move to final analysis
    if(DecisionInterimB1$decision=="Continue"){
        ResFinal <- AnalyzeData(d)
        DecisionFinalB1 <- Decision(ResFinal,  #results from AnalyzeData
                                    PlannedB1,  #results from CalcBoundaries
                                    k=2, #at which phase are we?
                                    analysis="final", #is it an interim or decision or final analysis
                                    Ik=c(ResInterim$Info,ResFinal$Info),  #all I_k (information) from first interim to final analysis (observed where possible)
                                    Id=NULL,  #expected or observed information at each decision analysis
                                    PositiveIsGood=FALSE, # whether positive effect is good (i.e. positive trial)
                                    Trace=FALSE, # whether to print some messages
                                    plot=FALSE)  #should the boundaries and results be plotted?
        ConclusionTrial.B1 <- DecisionFinalB1$decision
        critical.final.B1 <- DecisionFinalB1$details
        #--
        DecisionFinalB2 <- Decision(ResFinal,  #results from AnalyzeData
                                    PlannedB2,  #results from CalcBoundaries
                                    k=2, #at which phase are we?
                                    analysis="final", #is it an interim or decision or final analysis
                                    Ik=c(ResInterim$Info,ResFinal$Info),  #all I_k (information) from first interim to final analysis (observed where possible)
                                    Id=NULL,  #expected or observed information at each decision analysis
                                    PositiveIsGood=FALSE, # whether positive effect is good (i.e. positive trial)
                                    Trace=FALSE, # whether to print some messages
                                    plot=FALSE)  #should the boundaries and results be plotted?
        ConclusionTrial.B2 <- DecisionFinalB2$decision
        critical.final.B2 <- DecisionFinalB2$details
        #--                
        Z.final <- ResFinal$estimate/ResFinal$se
        est.final <- ResFinal$estimate
        se.final <- ResFinal$se
        Info.final <- ResFinal$Info
        PrCtInfo.final  <-  ResFinal$Info/PlannedB1$Imax
    }else{
        critical.final.B1 <- NA
        critical.final.B2 <- NA
        Z.final <- NA
        est.final <- NA
        se.final <- NA
        Info.final <- NA
        PrCtInfo.final  <- NA
    }    
    ## print(paste0("Conclusion=",ConclusionTrial.B1))    
    # }}}
    stopComp <- Sys.time()
    # {{{ Save results
    out <- c(ConclusionTrial.B1=ConclusionTrial.B1,
             ConclusionTrial.B2=ConclusionTrial.B2,
             seed=myseedi,             
             #---
             u.interim.B1 = u.interim.B1,
             l.interim.B1 = l.interim.B1,
             dec.interim.B1 = dec.interim.B1,
             c.decision.B1 = c.decision.B1,
             critical.final.B1 = critical.final.B1,
             #--
             u.interim.B2 = u.interim.B2,
             l.interim.B2 = l.interim.B2,
             dec.interim.B2 = dec.interim.B2,
             c.decision.B2 = c.decision.B2,
             critical.final.B2 = critical.final.B2,
             #--
             Z.interim = Z.interim,
             est.interim = est.interim,
             se.interim = se.interim,
             Info.interim = Info.interim,
             PrCtInfo.interim = PrCtInfo.interim,
             #--            
             Z.decision = Z.decision,
             est.decision = est.decision,
             se.decision = se.decision,
             Info.decision = Info.decision,
             PrCtInfo.decision = PrCtInfo.decision,
             #--            
             Z.final = Z.final,
             est.final = est.final,
             se.final = se.final,
             Info.final = Info.final,
             PrCtInfo.final = PrCtInfo.final,
             #----
             ## computation time
             CompTime=round(difftime(stopComp,startComp,units="secs"),3)
             )
    ## names(out) <- myColNames
    RES <- rbind(RES,out)
    # }}}
}
rownames(RES) <- NULL
save(RES,file=paste0(pathToSave,name,"-",i,".rda"))


# {{{ Quick look a the results
## RES
## summary(RES)
## head(RES[,order(colnames(RES))])
table(interim=RES[,"dec.interim.B1"],
      conclusion=RES[,"ConclusionTrial.B1"]
      )
table(interim=RES[,"dec.interim.B2"],
      conclusion=RES[,"ConclusionTrial.B2"]
      )
table(interim.B1=RES[,"dec.interim.B1"],
      interim.B2=RES[,"dec.interim.B2"]
      )
# }}}

#----------------------------------------------------------------------
### SimuMain.R ends here
