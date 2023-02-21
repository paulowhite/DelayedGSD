### BUILD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2022 (16:40) 
## Version: 
## Last-Updated: feb 20 2023 (14:27) 
##           By: Brice Ozenne
##     Update #: 71
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
path.Bresults <- file.path(path,"Results-built")
export <- TRUE

## * function used to collect the results from different files
loadRes <- function(path, tempo.file = FALSE, type = NULL,
                    export.attribute = NULL, trace = 2, space = "     "){
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
    if(n.file==0){
        if(trace>1){
            cat(space,"No file found in ",path,". \n\n",sep="")
        }
        return(NULL)
    }
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

    Urow <- unique(sapply(ls.out, NROW))
    if(length(Urow)==1){
        cat(space,length(ls.out)," files with ",Urow," lines \n\n",sep="")
    }else{
        cat(space,length(ls.out)," files, with from ",min(Urow)," to ",max(Urow)," lines \n\n",sep="")
    }
    return(out)
}

## * Process results
## ** Aggregate files and export 
dir.2stage <- grep("2stage",list.dirs(path = path.results), value = TRUE)
name.2stage <- stats::setNames(paste0("res2stage_",sapply(strsplit(dir.2stage, split = "2stage_"),"[",2)), dir.2stage)

for(iDir in dir.2stage){ ## iDir <- dir.2stage[1]
    iName <- name.2stage[iDir]
    cat(which(iDir == dir.2stage),") read \'",iDir,"\' \n",
        "     save in \'",iName,"\' \n", sep = "")
    assign(x = iName, value = loadRes(iDir))
    if(export){
        saveRDS(eval(parse(text=iName)), file = file.path(path.results,paste0(iName,".rds") ))
    }
}


## dt <- loadRes("x:/DelayedGSD/Results/2stage_missing_binding_ar10_power", tempo.file = TRUE)
## dt[, .N, by = "file"]
## length(unique(dt$file))
## dt[file=="sim-2stage_missing_binding_ar10_power-1(tempo)_100.rds"]


## length(unique(res2stage_missing_binding_ar10_power$file)) ## 99
## length(unique(res2stage_missing_binding_ar10_typeI$file))        ## 99
## length(unique(res2stage_missing_binding_ar5_power$file))         
## length(unique(res2stage_missing_binding_ar5_typeI$file))      
## length(unique(res2stage_missing_fixC_binding_ar10_power$file)) ## 98
## length(unique(res2stage_missing_fixC_binding_ar10_typeI$file))   ## 98
## length(unique(res2stage_missing_fixC_binding_ar5_power$file))    
## length(unique(res2stage_missing_fixC_binding_ar5_typeI$file))    
## length(unique(res2stage_missing_fixC_nonbinding_ar10_power$file)) ## 99
## length(unique(res2stage_missing_fixC_nonbinding_ar10_typeI$file)) ## 99
## length(unique(res2stage_missing_fixC_nonbinding_ar5_power$file)) 
## length(unique(res2stage_missing_fixC_nonbinding_ar5_typeI$file)) 
## length(unique(res2stage_missing_nonbinding_ar10_power$file))     
## length(unique(res2stage_missing_nonbinding_ar10_typeI$file))     
## length(unique(res2stage_missing_nonbinding_ar5_power$file))      
## length(unique(res2stage_missing_nonbinding_ar5_typeI$file))      
## length(unique(res2stage_nomissing_binding_ar5_power$file))       
## length(unique(res2stage_nomissing_binding_ar5_typeI$file))       

## ** Aggregate scenario and export 
legend.2stage <- data.frame(name = name.2stage,
                            scenario = 1:length(name.2stage),
                            missing = grepl("nomissing", name.2stage)== FALSE,
                            binding = grepl("nonbinding", name.2stage)== FALSE,
                            fixC = grepl("fixC", name.2stage),
                            ar = gsub("ar","",sapply(strsplit(name.2stage, split = "_"), function(iVec){tail(iVec,2)[1]})),
                            hypo = sapply(strsplit(name.2stage, split = "_"), function(iVec){tail(iVec,1)})
                            )

res2stage <- do.call(rbind,lapply(name.2stage, function(iName){ ## iName <- name.Bresults[1]
    data.table(legend.2stage[name.2stage==iName,-1], eval(parse(text = iName)))
}))
res2stage$computation.time <- NULL
res2stage$file <- NULL
res2stage$sigma <- NULL
res2stage$ar <- as.numeric(res2stage$ar)

if(export){
    saveRDS(res2stage, file = file.path(path.Bresults,"res2stage.rds"))
}

## res2stage <- readRDS(file = file.path(path.Bresults,"res2stage.rds"))
## unique(res2stage[,.N, by = c("scenario","method","type")]$N)
## [1] 10000


## * Analyse results (COBA)
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
