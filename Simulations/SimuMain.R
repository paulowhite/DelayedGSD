### simumain.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar  5 2021 (10:56) 
## Version: 
## Last-Updated: feb 18 2022 (11:47) 
##           By: Brice Ozenne
##     Update #: 478
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
## * parameters
name <- "ScenarioName" # To save the results
                                        #----
NMC <- 10 # number of sequential simulations to run in parallel. Eg. with 250, then we can run 40 scripts in paralell to get N=10,000 runs in total.
                                        #---
myseed <- 140786598
                                        #--- to plan the trial ----
kMax <- 2  #max number of analyses (including final)
alpha <- 0.025  #type I error (one sided)
beta <- 0.2  #type II error
informationRates <- c(0.5,1)  #planned  information rates
rho_alpha <- 2  # rho parameter for alpha error spending function
rho_beta <- 2  # rho parameter for beta error spending function
## deltaPower <- 0.75 # just to try another value when Id > Imax
Id <- 0.55  #(expected) information rate at each decision analysis
binding <- TRUE
                                        #
                                        #---- to generate data -----------
                                        #
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

# {{{ set path to load and save and othe machine specific variables
## * path
if(system("whoami",intern=TRUE)=="paul"){  
    pathToLoad <- "~/research/SeqDesignDelayed/DelayedGSD/R/"
    pathToSave <- "~/research/SeqDesignDelayed/DelayedGSD/Simulations/output/"
    name <- paste0(name,"test")
    i <- 1
}else if(system("whoami",intern=TRUE) %in% c("bozenne","unicph\\hpl802")){  
    pathToLoad <- "~/Documents/GitHub/DelayedGSD/R"
    i <- 1
}else if(system("whoami",intern=TRUE)=="tfq625"){
    pathToLoad <- "~/to-sync-server-biostat/TOCREATE"
    pathToSave <- "~/to-sync-server-biostat/TOCREATE"
    i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}
# }}}


## * libraries and functions
library(DelayedGSD)
## sourceDir <- function(path, trace = TRUE, ...) {
##     for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
##         if(trace) cat(nm,":")
##         source(file.path(path, nm), ...)
##         if(trace) cat("\n")
##     }
## }
## sourceDir(pathToLoad)
    

## * Compute inflation factor and sample size
plannedB <- vector(mode = "list", length = 3)
for(iMeth in 1:3){ ## iMeth <- 1
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
                                        delta=tail(delta,1))
    ## summary(plannedB[[1]])
    ## coef(plannedB[[iMeth]], type = "information")
}
inflationFactor <- sapply(plannedB,function(iP){iP$planned$InflationFactor})
nGSD <- ceiling(n*inflationFactor)
##  plot(plannedB[[1]])

#
# --- just to check---
## n
## InfoFixed <- ((qnorm(1-beta)-qnorm(alpha))/(deltaPower))^2
## InfoFixed
## n/(4*(allsd[3])^2)
#---
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


for(j in allj){ ## j <- 5 ## 5

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
    ## ** simulate
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
    thet <- d$t3[ceiling(n*PropForInterim)]
    di <- SelectData(d,t=thet)
    ## ddi <- FormatAsCase(di) # needed ????
    ## head(d[d$id==52,])
                                        # }}}
    ## {{{ analyze data at at interim
    ## ** interim
    lmmI <- analyzeData(di, ddf = "nlme", data.decision = sum(d$t1 <= thet + theDelta.t*TimeFactor), getinfo = TRUE, trace = TRUE)
    ## lmmI <- analyzeData(di, ddf = "nlme", getinfo = TRUE, trace = TRUE)

    currentGSD <- vector(mode = "list", length = 3)
    out.interim <- vector(mode = "list", length = 3)
    for(iMeth in 1:3){ ## iMeth <- 1

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
                                            reason = iDecision.interim["reason.interim","stage 1"])
    }
 ## currentGSD[[1]]

    ## ** decision
    dDecision <- d[which(d$t1 <= thet + theDelta.t*TimeFactor),]
    lmmD <- analyzeData(dDecision, ddf = "nlme", getinfo = TRUE, trace = TRUE)
    
    out.decision <- vector(mode = "list", length = 3)
    for(iMeth in 1:3){ ## iMeth <- 1
          
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
                                                info = iInfo.decision[1,"Interim"],
                                                infoPC = iInfo.decision[1,"Interim.pc"],
                                                ck = iBoundary.decision[1,"Cbound"],
                                                decision = unname(iDecision.decision["decision","stage 2"])
                                                )

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
                                                decision = NA)
        }
    }
                                        # }}}
                                        # {{{ Analyze data at decision

    ## ** finale
    dFinal <- d
    lmmF <- analyzeData(dFinal, ddf = "nlme", getinfo = TRUE, trace = TRUE)

    out.final <- c()
    for(iMeth in 1:3){ ## iMeth <- 1
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
                                             decision = NA)

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
                                             info = iInfo.final[1,"Interim"],
                                             infoPC = iInfo.final[1,"Interim.pc"],
                                             ck = iBoundary.final[1,"Cbound"],
                                             decision = unname(coef(currentGSD[[iMeth]], type = "decision")["decision","stage 2"])
                                             )
        }
    }
                                        # }}}

    stopComp <- Sys.time()
                                        # {{{ Save results

    outMerge <- do.call(rbind,lapply(1:3, function(iMeth){
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
        time.interim = thet,
        seed=myseedi,             
        nX1.interim = sum(!is.na(di$X1)),
        nX2.interim = sum(!is.na(di$X2)),
        nX3.interim = sum(!is.na(di$X3)),
        ## computation time
        computation.time=as.double(round(difftime(stopComp,startComp,units="secs"),3))
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
