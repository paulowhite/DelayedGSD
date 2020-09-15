### main.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 24 2020 (14:37) 
## Version: 
## Last-Updated: Sep 15 2020 (11:33) 
##           By: Paul Blanche
##     Update #: 169
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
source("Rfunctions/GenData.R")
source("Rfunctions/PlotProgress.R")
source("Rfunctions/AnalyzeData.R")
source("Rfunctions/blockrand.R")
source("Rfunctions/SelectData.R")


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



#----------------------------------------------------------------------
### main.R ends here
