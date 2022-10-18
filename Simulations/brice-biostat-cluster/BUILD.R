### BUILD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2022 (16:40) 
## Version: 
## Last-Updated: okt 18 2022 (13:58) 
##           By: Brice Ozenne
##     Update #: 16
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

## * function used to collect the results from different files
loadRes <- function(path, tempo.file = FALSE, type = NULL,
                    export.attribute = NULL, trace = TRUE){
    all.files <- list.files(path)
    file.tempo <- grep("(tempo)",all.files,value = TRUE)
    file.final <- setdiff(all.files, file.tempo)

    if(tempo.file){
        file.read <- file.tempo
    }else{
        file.read <- file.final
    }
    if(!is.null(type)){
        file.read <- grep(pattern=type,x=file.read,value=TRUE)
    }

    n.file <- length(file.read)

    myApply <- switch(as.character(as.logical(trace)),
                      "TRUE" = pbapply::pblapply,
                      "FALSE" = lapply)

    ls.out <- do.call(myApply, args = list(X = 1:n.file, FUN = function(iFile){
        if(grepl("\\.rds$",file.read[iFile])){
            iRead <- try(readRDS(file = file.path(path,file.read[iFile])))
        }else if(grepl("\\.rda$", file.read[iFile])){
            iRead <- try(load(file = file.path(path,file.read[iFile])))
            if(!inherits(iRead,"try-error")){
                iRead <- eval(parse(text = iRead))
            }
        }else{
            return(NULL)
        }
        if(inherits(iRead,"try-error")){
            return(NULL)
        }else{
            iOut <- cbind(data.table::as.data.table(iRead),
                          file = file.read[iFile])
            return(iOut)
        }
    }))
    out <- do.call(rbind, ls.out)
    return(out)
}

## * Collect files
dt.missing_power <- loadRes(file.path(path.results,"missing_power"))
dt.missing_typeI <- loadRes(file.path(path.results,"missing_typeI"))
dt.noFixC_power <- loadRes(file.path(path.results,"noFixC_power"))
dt.noFixC_typeI <- loadRes(file.path(path.results,"noFixC_typeI"))
dt.nomissing_power <- loadRes(file.path(path.results,"nomissing_power"))
dt.nomissing_typeI <- loadRes(file.path(path.results,"nomissing_typeI"))
dt.nonbinding_power <- loadRes(file.path(path.results,"nonbinding_power"))
dt.nonbinding_typeI <- loadRes(file.path(path.results,"nonbinding_typeI"))

## * Export
saveRDS(dt.missing_power, file = file.path(path.results,"res-simuCOBA-missing_power.rds") )
saveRDS(dt.missing_typeI, file = file.path(path.results,"res-simuCOBA-missing_typeI.rds") )
saveRDS(dt.noFixC_power, file = file.path(path.results,"res-simuCOBA-noFixC_power.rds") )
saveRDS(dt.noFixC_typeI, file = file.path(path.results,"res-simuCOBA-noFixC_typeI.rds") )
saveRDS(dt.nomissing_power, file = file.path(path.results,"res-simuCOBA-nomissing_power.rds") )
saveRDS(dt.nomissing_typeI, file = file.path(path.results,"res-simuCOBA-nomissing_typeI.rds") )
saveRDS(dt.nonbinding_power, file = file.path(path.results,"res-simuCOBA-nonbinding_power.rds") )
saveRDS(dt.nonbinding_typeI, file = file.path(path.results,"res-simuCOBA-nonbinding_typeI.rds") )

## * Inspect results
## dt.missing_power <- readRDS(file = file.path(path.results,"res-simuCOBA-missing.rds") )

## table(table(dt.sim$seed))
dt.missing_power[,.(N = .N, results = "efficacy" %in% na.omit(decision)),by = c("method","seed")][,.(power=mean(results)),by="method"]
##    method  power
## 1:      1 0.8274
## 2:      2 0.8273
## 3:      3 0.8266

dcast(dt.sim[type!="final",.(N=.N,infoPC=mean(infoPC)),by=c("method","type")],
      N+method~type, value.var = "infoPC")
##        N method  decision   interim
## 1: 10000      1 0.5630846 0.5199598
## 2: 10000      2 0.5630846 0.5199598
## 3: 10000      3 0.5641234 0.5198136

dt.sim[p.value_MUE<0.025 & decision=="futility",.(method,stage,type,seed)]
##    method stage  type      seed
## 1:      3     2 final 65836753
dt.sim[p.value_MUE>0.025 & decision=="efficacy"]

##----------------------------------------------------------------------
### BUILD.R ends here
