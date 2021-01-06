### main.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 24 2020 (14:37) 
## Version: 
## Last-Updated: Jan  6 2021 (14:29) 
##           By: Paul Blanche
##     Update #: 190
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

rm(list=ls())

if(system("whoami",intern=TRUE)=="paul"){  
    setwd("~/research/SeqDesignDelayed/DelayedGSD/")
}
sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}
sourceDir("Rfunctions/")


# generate data
res <- GenData(n=100,N.fw=4,ar=10)
d <-res$d
## res$cor
## head(d)
## tail(d)
#
#


#----- select when half of the subjects have one follow-up measuement---
# which depends on accrual rate (ar) and time to process data (Delta.t)
theDelta.t <- 0.500001
thear <- 10
thet <- d$t2[ceiling(nrow(d)/2) + ceiling(thear*theDelta.t)]

# plot progress in including patients and collected data
PlotProgress(d,at=thet)

# Create data available for nalaysis at that time
di <- SelectData(d,t=thet)
di
# Analyze the data
#
Res <- AnalyzeData(di)
Res$estimate
Res$se
summary(Res$fit)

#--- plan boundaries ---
PlannedB <- CalcBoundaries(kMax=2,  #max number of analyses (including final)
                           sided=1,  #one or two-sided
                           alpha=0.025,  #type I error
                           beta=0.2,  #type II error
                           informationRates=c(0.5,1),  #planned or observed information rates
                           gammaA=2,  #rho parameter for alpha error spending function
                           gammaB=2,  #rho parameter for beta error spending function
                           method=1,  #use method 1 or 2 from paper H&J
                           cNotBelowFixedc=TRUE, # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                           delta=1.5,  #effect that the study is powered for
                           Id=0.55)    #(expected) information ratio at each decision analysis

#-----
# Here is a weird example which seems to indicate a BUG !!!!!!!!!!!!!!!!!!!
#-----
## PlannedB <- CalcBoundaries(kMax=3,  #max number of analyses (including final)
                           ## sided=1,  #one or two-sided
                           ## alpha=0.025,  #type I error
                           ## beta=0.2,  #type II error
                           ## informationRates=c(0.5,0.75,1),  #planned or observed information rates
                           ## gammaA=1,  #rho parameter for alpha error spending function
                           ## gammaB=1,  #rho parameter for beta error spending function
                           ## method=1,  #use method 1 or 2 from paper H&J
                           ## cNotBelowFixedc=TRUE, # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                           ## delta=1.5,  #effect that the study is powered for
                           ## Id=c(0.55,0.85))    #(expected) information ratio at each decision analysis

PlannedB

par(mfrow=c(1,3))
PlotBoundaries(PlannedB,type="Z",Itype="rate")
PlotBoundaries(PlannedB,type="P")
PlotBoundaries(PlannedB,type="E",Itype="abs")


PlannedB1F <- CalcBoundaries(method=1,  #use method 1 or 2 from paper H&J
                             cNotBelowFixedc=FALSE # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                             )
PlannedB2F <- CalcBoundaries(method=2,  #use method 1 or 2 from paper H&J
                             cNotBelowFixedc=FALSE # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                             )    
PlannedB1T <- CalcBoundaries(method=1,  #use method 1 or 2 from paper H&J
                             cNotBelowFixedc=TRUE # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                             )
PlannedB2T <- CalcBoundaries(method=2,  #use method 1 or 2 from paper H&J
                             cNotBelowFixedc=TRUE # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                             )    

par(mfrow=c(2,2))
PlotBoundaries(PlannedB1F,type="Z",Itype="rate")
title("Method 1, c value at decision allowed too low")
PlotBoundaries(PlannedB2F,type="Z",Itype="rate")
title("Method 2, c value at decision allowed too low")
PlotBoundaries(PlannedB1T,type="Z",Itype="rate")
title("Method 1, c value at decision NOT allowed too low")
PlotBoundaries(PlannedB2T,type="Z",Itype="rate")
title("Method 2, c value at decision NOT allowed too low")


#----------------------------------------------------------------------
### main.R ends here

