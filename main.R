### main.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 24 2020 (14:37) 
## Version: 
## Last-Updated: feb  3 2021 (14:06) 
##           By: Brice Ozenne
##     Update #: 173
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
}else if(system("whoami",intern=TRUE)=="brice"){  
    setwd("~/Documents/GitHub/DelayedGSD/")
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

Res <- AnalyzeData(di)
Res <- AnalyzeData(di, ddf = "satterthwaite")


ResInfo <- getInformation(Res$fit, data = Res$d.long, name.coef = "Z1", type = "estimation", details = TRUE)
as.double(ResInfo) - 1/vcov(Res$fit)["Z1","Z1"]
attr(ResInfo,"details")$decision$information
attr(ResInfo,"details")$interim$information
attr(ResInfo,"details")$interim.cc$information

getInformation(Res$fit, data = Res$d.long, name.coef = "Z1", type = "prediction", method.prediction = "inflation")




#----------------------------------------------------------------------
### main.R ends here
