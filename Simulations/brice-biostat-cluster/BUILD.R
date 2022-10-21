### BUILD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2022 (16:40) 
## Version: 
## Last-Updated: okt 21 2022 (10:29) 
##           By: Brice Ozenne
##     Update #: 26
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

## ** power
dt.missing_power[,.(N = .N, results = "efficacy" %in% na.omit(decision)),by = c("method","seed")][,.(power=mean(results)),by="method"]
##    method  power
## 1:      1 0.8060
## 2:      2 0.8045
## 3:      3 0.8021

dt.noFixC_power[,.(N = .N, results = "efficacy" %in% na.omit(decision)),by = c("method","seed")][,.(power=mean(results)),by="method"]
##    method  power
## 1:      1 0.8000
## 2:      2 0.8046
## 3:      3 0.8008

dt.nomissing_power[,.(N = .N, results = "efficacy" %in% na.omit(decision)),by = c("method","seed")][,.(power=mean(results)),by="method"]
##    method  power
## 1:      1 0.8031
## 2:      2 0.8028
## 3:      3 0.7993

dt.nonbinding_power[,.(N = .N, results = "efficacy" %in% na.omit(decision)),by = c("method","seed")][,.(power=mean(results)),by="method"]
##    method  power
## 1:      1 0.8065
## 2:      2 0.8064
## 3:      3 0.8046

## ** type 1 error
dt.missing_typeI[,.(N = .N, results = "efficacy" %in% na.omit(decision)),by = c("method","seed")][,.(power=mean(results)),by="method"]
##    method  power
## 1:      1 0.0242
## 2:      2 0.0241
## 3:      3 0.0240

dt.noFixC_typeI[,.(N = .N, results = "efficacy" %in% na.omit(decision)),by = c("method","seed")][,.(power=mean(results)),by="method"]
##    method  power
## 1:      1 0.0242
## 2:      2 0.0239
## 3:      3 0.0250

dt.nomissing_typeI[,.(N = .N, results = "efficacy" %in% na.omit(decision)),by = c("method","seed")][,.(power=mean(results)),by="method"]
##    method  power
## 1:      1 0.0246
## 2:      2 0.0244
## 3:      3 0.0245

dt.nonbinding_typeI[,.(N = .N, results = "efficacy" %in% na.omit(decision)),by = c("method","seed")][,.(power=mean(results)),by="method"]
##    method  power
## 1:      1 0.0239
## 2:      2 0.0238
## 3:      3 0.0246

## ** special cases
table(dt.missing_power$reason)
## decreasing information               efficacy               futility    no boundary crossed 
##                      3                  10825                   1793                  17382 
table(dt.noFixC_power$reason)
## efficacy            futility no boundary crossed 
##    10680                1594               17726 
table(dt.nomissing_power$reason)
## efficacy            futility no boundary crossed 
##    10289                1551               18160 
table(dt.nonbinding_power$reason)
## efficacy            futility no boundary crossed 
##    11139                1679               17182 

table(dt.missing_power$reason)
## decreasing information               efficacy               futility    no boundary crossed 
##                      3                  10825                   1793                  17382 
table(dt.noFixC_power$reason)
## efficacy            futility no boundary crossed 
##    10680                1594               17726 
table(dt.nomissing_power$reason)
## efficacy            futility no boundary crossed 
##    10289                1551               18160 
table(dt.nonbinding_power$reason)
## efficacy            futility no boundary crossed 
##    11139                1679               17182 

## ** reversal
dtRev.missing_power <- dt.missing_power[,.(reversal = (decision[1] == "stop") * ((reason[1] == "efficacy" & decision[2] == "futility")+(reason[1] == "futility" & decision[2] == "efficacy"))) ,
                                        by = c("method","seed")]
dtRev.missing_power$method


dcast(dt.sim[type!="final",.(N=.N,infoPC=mean(infoPC)),by=c("method","type")],
      N+method~type, value.var = "infoPC")
##        N method  decision   interim
## 1: 10000      1 0.5630846 0.5199598
## 2: 10000      2 0.5630846 0.5199598
## 3: 10000      3 0.5641234 0.5198136

## ** frequency mismatch p-value / boundaries
dt.sim[p.value_MUE<0.025 & decision=="futility",.(method,stage,type,seed)]
##    method stage  type      seed
## 1:      3     2 final 65836753
dt.sim[p.value_MUE>0.025 & decision=="efficacy"]

## ** bias (ML, MUE)

## ** Process as COBA
res <- dt.noFixC_power
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
result



discrep <- res[res$p.value_MUE>0.025 & !is.na(res$p.value_MUE),]
discrep <- discrep[!is.na(discrep$final.efficacy),]
dim(discrep)

discrep <- res[res$p.value_MUE<0.025 & !is.na(res$p.value_MUE),]
discrep <- discrep[is.na(discrep$final.efficacy),]
dim(discrep)
table(discrep$method)
## length(unique(paste(discrep$seed,discrep$method)))

## refinement: frequency stratified by method

res_temp[res_temp$seed==res_temp$seed[2],]

index.NNA <- which(!is.na(res_temp$estimate_MUE))
myestimate_MUE <- res_temp[index.NNA,estimate_MUE]
myestimate_ML <- res_temp[index.NNA,estimate_ML]
dt.gg <- res_temp[index.NNA,.(estimate_MUE,estimate_ML,stage)]
dt.gg$stage <- as.factor(dt.gg$stage)
library(ggplot2)
ggplot(dt.gg, aes(x = estimate_ML, group = stage, color = stage, fill = stage)) + geom_histogram() + facet_wrap(~stage)


mean(myestimate_MUE)
mean(myestimate_ML)
hist(myestimate_MUE, breaks = 20)
hist(myestimate_ML, breaks = 20)
hist(myestimate_MUE-myestimate_ML, breaks = 20)

##----------------------------------------------------------------------
### BUILD.R ends here
