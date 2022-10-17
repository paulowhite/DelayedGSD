### BUILD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2022 (16:40) 
## Version: 
## Last-Updated: okt 13 2022 (11:23) 
##           By: Brice Ozenne
##     Update #: 7
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(data.table)
library(ggplot2)

path <- "x:/DelayedGSD"
path.results <- file.path(path,"Results")

## * Collect files
vec.file <- list.files(file.path(path.results."SimuCOBA-light"))

ls.file <- lapply(vec.file, function(iFile){
    load(file.path(path.results."SimuCOBA-light",iFile))
    return(RES)
})
dt.sim <- as.data.table(do.call(rbind,ls.file))
## export
saveRDS(dt.sim, file = file.path(path.results,"res-simuCOBA-light.rds") )

## * Collect files
## dt.sim <- readRDS(file = file.path(path.results,"res-simuCOBA-light.rds") )

## table(table(dt.sim$seed))
dt.sim[,.(N = .N, results = "efficacy" %in% na.omit(decision)),by = c("method","seed")][,.(power=mean(results)),by="method"]
##    method  power
## 1:      1 0.8009
## 2:      2 0.8009
## 3:      3 0.7979

dcast(dt.sim[type!="final",.(N=.N,infoPC=mean(infoPC)),by=c("method","type")],
      N+method~type, value.var = "infoPC")
##        N method  decision   interim
## 1: 10000      1 0.5630846 0.5199598
## 2: 10000      2 0.5630846 0.5199598
## 3: 10000      3 0.5641234 0.5198136

dt.sim[p.value_MUE<0.025 & decision=="futility",.(method,stage,type,seed)]
##    method stage  type      seed
## 1:      3     2 final  65836753
## 2:      3     2 final 522862069
## 3:      3     2 final 632796957
## 4:      3     2 final 502705982
## 5:      3     2 final 545711051
## 6:      3     2 final 866677443
## 7:      3     2 final 524781539
## 8:      3     2 final 782550339
## 9:      3     2 final 826717444
dt.sim[p.value_MUE>0.025 & decision=="efficacy"]

##----------------------------------------------------------------------
### BUILD.R ends here
