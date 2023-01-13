### BUILD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2022 (16:40) 
## Version: 
## Last-Updated: jan 13 2023 (13:39) 
##           By: Brice Ozenne
##     Update #: 64
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

## * Process results
## ** Aggregate files and export 
dir.2stage <- grep("2stage",list.dirs(path = path.results), value = TRUE)
name.2stage <- stats::setNames(paste0("res2stage_",sapply(strsplit(dir.2stage, split = "2stage_"),"[",2)), dir.2stage)

for(iDir in dir.2stage){ ## iDir <- dir.2stage[1]
    iName <- name.2stage[iDir]
    cat(which(iDir == dir.2stage),") read \'",iDir,"\' and save in \'",iName,"\' \n", sep = "")
    assign(x = iName, value = loadRes(iDir))
    if(export){
        saveRDS(eval(parse(text=iName)), file = file.path(path.results,paste0(iName,".rds") ))
    }
}


## dt <- loadRes("x:/DelayedGSD/Results/2stage_missing_binding_ar10_power", tempo.file = TRUE)
## dt[, .N, by = "file"]
## length(unique(dt$file))
## dt[file=="sim-2stage_missing_binding_ar10_power-1(tempo)_100.rds"]

dt <- loadRes("x:/DelayedGSD/Results/2stage_missing_fixC_binding_ar10_typeI", tempo.file = TRUE)
table(dt[,.N,by = "file"]$N)
## 55:  sim-2stage_missing_fixC_binding_ar10_typeI-59(tempo)_100.rds 630


length(unique(res2stage_missing_binding_ar10_power$file)) ## 99
length(unique(res2stage_missing_binding_ar10_typeI$file))        ## 99
## length(unique(res2stage_missing_binding_ar5_power$file))         
## length(unique(res2stage_missing_binding_ar5_typeI$file))      
length(unique(res2stage_missing_fixC_binding_ar10_power$file)) ## 98
length(unique(res2stage_missing_fixC_binding_ar10_typeI$file))   ## 98
## length(unique(res2stage_missing_fixC_binding_ar5_power$file))    
## length(unique(res2stage_missing_fixC_binding_ar5_typeI$file))    
length(unique(res2stage_missing_fixC_nonbinding_ar10_power$file)) ## 99
length(unique(res2stage_missing_fixC_nonbinding_ar10_typeI$file)) ## 99
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

## * Analyse results (data.table)
## ** rejection rate
## for each run, create a binary indicator for rejection for efficacy
res2stage.rejection <- res2stage[missing==TRUE,.(N = .N, rejection = "efficacy" %in% na.omit(decision)),by = c("method","seed","binding","fixC","ar","hypo")]
## all(res2stage.rejection$N==3)
## stage 1 (interim), stage 1 (decision), stage 2

## average over runs and method within scenario
res2stageS.rejection <- res2stage.rejection[,.(N = .N, rejectionRate = 100*mean(rejection)),by=c("method","binding","fixC","ar","hypo")]

## power by by method (columns) and scenario (rows)
dcast(res2stageS.rejection[hypo=="power"], binding + fixC + ar ~ method, value.var = "rejectionRate")
##    binding  fixC ar     1     2     3
## 1:   FALSE FALSE  5 80.65 80.64 80.46
## 2:   FALSE FALSE 10 80.67 80.60 80.65
## 3:   FALSE  TRUE  5 80.65 80.64 80.46
## 4:   FALSE  TRUE 10 80.67 80.60 80.65
## 5:    TRUE FALSE  5 80.60 80.45 80.21
## 6:    TRUE FALSE 10 81.00 80.79 80.45
## 7:    TRUE  TRUE  5 80.60 80.45 80.21
## 8:    TRUE  TRUE 10 81.00 80.79 80.45

## type 1 error by method (columns) and scenario (rows)
dcast(res2stageS.rejection[hypo=="typeI"], binding + fixC + ar ~ method, value.var = "rejectionRate")
##    binding  fixC ar    1    2    3
## 1:   FALSE FALSE  5 2.39 2.38 2.46
## 2:   FALSE FALSE 10 2.35 2.33 2.47
## 3:   FALSE  TRUE  5 2.39 2.38 2.46
## 4:   FALSE  TRUE 10 2.35 2.33 2.47
## 5:    TRUE FALSE  5 2.42 2.41 2.40
## 6:    TRUE FALSE 10 2.46 2.53 2.40
## 7:    TRUE  TRUE  5 2.42 2.41 2.40
## 8:    TRUE  TRUE 10 2.46 2.53 2.40

## ** Conclusion of the trial
res2stageS.final <- res2stage[missing == TRUE & decision %in% c("futility","efficacy"),
                             .(N = .N,
                               interim.efficacy = 100*mean((stage == 1)*(decision == "efficacy")),
                               interim.futility = 100*mean((stage == 1)*(decision == "futility")),
                               final.efficacy = 100*mean((stage == 2)*(decision == "efficacy")),
                               final.futility = 100*mean((stage == 2)*(decision == "futility"))),
                             by = c("method","binding","fixC","ar","hypo")]
dcast(res2stageS.final[method==1], hypo + binding + fixC + ar ~ method, value.var = c("interim.efficacy","interim.futility","final.efficacy","final.futility"))
##      hypo binding  fixC ar interim.efficacy_1 interim.futility_1 final.efficacy_1 final.futility_1
##  1: power   FALSE FALSE  5              36.75               5.70            43.90            13.65
##  2: power   FALSE FALSE 10              38.32               5.87            42.35            13.46
##  3: power   FALSE  TRUE  5              36.75               5.70            43.90            13.65
##  4: power   FALSE  TRUE 10              38.32               5.87            42.35            13.46
##  5: power    TRUE FALSE  5              35.60               6.02            45.00            13.38
##  6: power    TRUE FALSE 10              37.82               6.05            43.18            12.95
##  7: power    TRUE  TRUE  5              35.60               6.02            45.00            13.38
##  8: power    TRUE  TRUE 10              37.82               6.05            43.18            12.95
##  9: typeI   FALSE FALSE  5               0.68              69.19             1.71            28.42
## 10: typeI   FALSE FALSE 10               0.75              71.29             1.60            26.36
## 11: typeI   FALSE  TRUE  5               0.68              69.19             1.71            28.42
## 12: typeI   FALSE  TRUE 10               0.75              71.29             1.60            26.36
## 13: typeI    TRUE FALSE  5               0.68              69.21             1.74            28.37
## 14: typeI    TRUE FALSE 10               0.79              70.85             1.67            26.69
## 15: typeI    TRUE  TRUE  5               0.68              69.21             1.74            28.37
## 16: typeI    TRUE  TRUE 10               0.79              70.85             1.67            26.69

dcast(res2stageS.final[method==2], hypo + binding + fixC + ar ~ method, value.var = c("interim.efficacy","interim.futility","final.efficacy","final.futility"))
##      hypo binding  fixC ar interim.efficacy_2 interim.futility_2 final.efficacy_2 final.futility_2
##  1: power   FALSE FALSE  5              36.78               5.72            43.86            13.64
##  2: power   FALSE FALSE 10              38.33               6.11            42.27            13.29
##  3: power   FALSE  TRUE  5              36.78               5.72            43.86            13.64
##  4: power   FALSE  TRUE 10              38.33               6.11            42.27            13.29
##  5: power    TRUE FALSE  5              35.55               6.10            44.90            13.45
##  6: power    TRUE FALSE 10              37.66               6.22            43.13            12.99
##  7: power    TRUE  TRUE  5              35.55               6.10            44.90            13.45
##  8: power    TRUE  TRUE 10              37.66               6.22            43.13            12.99
##  9: typeI   FALSE FALSE  5               0.67              69.25             1.71            28.37
## 10: typeI   FALSE FALSE 10               0.74              71.84             1.59            25.83
## 11: typeI   FALSE  TRUE  5               0.67              69.25             1.71            28.37
## 12: typeI   FALSE  TRUE 10               0.74              71.84             1.59            25.83
## 13: typeI    TRUE FALSE  5               0.67              69.05             1.74            28.54
## 14: typeI    TRUE FALSE 10               0.85              71.18             1.68            26.29
## 15: typeI    TRUE  TRUE  5               0.67              69.05             1.74            28.54
## 16: typeI    TRUE  TRUE 10               0.85              71.18             1.68            26.29

dcast(res2stageS.final[method==3], hypo + binding + fixC + ar ~ method, value.var = c("interim.efficacy","interim.futility","final.efficacy","final.futility"))
##      hypo binding  fixC ar interim.efficacy_3 interim.futility_3 final.efficacy_3 final.futility_3
##  1: power   FALSE FALSE  5              37.37               5.86            43.09            13.68
##  2: power   FALSE FALSE 10              41.47               6.05            39.18            13.30
##  3: power   FALSE  TRUE  5              37.37               5.86            43.09            13.68
##  4: power   FALSE  TRUE 10              41.47               6.05            39.18            13.30
##  5: power    TRUE FALSE  5              36.49               6.42            43.72            13.37
##  6: power    TRUE FALSE 10              40.44               6.54            40.01            13.01
##  7: power    TRUE  TRUE  5              36.49               6.42            43.72            13.37
##  8: power    TRUE  TRUE 10              40.44               6.54            40.01            13.01
##  9: typeI   FALSE FALSE  5               0.75              68.45             1.71            29.09
## 10: typeI   FALSE FALSE 10               0.82              69.00             1.65            28.53
## 11: typeI   FALSE  TRUE  5               0.75              68.45             1.71            29.09
## 12: typeI   FALSE  TRUE 10               0.82              69.00             1.65            28.53
## 13: typeI    TRUE FALSE  5               0.68              68.37             1.72            29.23
## 14: typeI    TRUE FALSE 10               0.74              68.77             1.66            28.83
## 15: typeI    TRUE  TRUE  5               0.68              68.37             1.72            29.23
## 16: typeI    TRUE  TRUE 10               0.74              68.77             1.66            28.83

## ** bias
true_eff <- 0.6
## for each run, error made by each estimator
res2stage$truth <- c(0,true_eff)[(res2stage$hypo=="power")+1]
res2stage.bias <- res2stage[missing == TRUE & decision %in% c("futility","efficacy"),
                            .(N = .N, bias_MLE = estimate_ML-truth, bias_MUE = (estimate_MUE>truth) - 0.5),
                            by = c("method","seed","binding","fixC","ar","hypo")]
res2stageS.bias <- res2stage.bias[,.(N = .N, bias_MLE = mean(bias_MLE, na.rm = TRUE), bias_MUE = mean(bias_MUE, na.rm = TRUE)),by=c("method","binding","fixC","ar","hypo")]

dcast(res2stageS.bias[hypo=="typeI"], binding + fixC + ar ~ method, value.var = c("bias_MLE","bias_MUE"))
##    binding  fixC ar  bias_MLE_1  bias_MLE_2  bias_MLE_3 bias_MUE_1 bias_MUE_2 bias_MUE_3
## 1:   FALSE FALSE  5 -0.03190984 -0.03193420 -0.03210110    -0.0544    -0.0543    -0.0571
## 2:   FALSE FALSE 10 -0.02007496 -0.01979455 -0.02024580    -0.0458    -0.0446    -0.0524
## 3:   FALSE  TRUE  5 -0.03190984 -0.03193420 -0.03210110    -0.0544    -0.0543    -0.0571
## 4:   FALSE  TRUE 10 -0.02007496 -0.01979455 -0.02024580    -0.0458    -0.0446    -0.0524
## 5:    TRUE FALSE  5 -0.03041945 -0.03082210 -0.03057683     0.0000    -0.0002     0.0001
## 6:    TRUE FALSE 10 -0.01841568 -0.01842973 -0.01850869     0.0002    -0.0013     0.0001
## 7:    TRUE  TRUE  5 -0.03041945 -0.03082210 -0.03057683     0.0000    -0.0002     0.0001
## 8:    TRUE  TRUE 10 -0.01841568 -0.01842973 -0.01850869     0.0002    -0.0013     0.0001

dcast(res2stageS.bias[hypo=="power"], binding + fixC + ar ~ method, value.var = c("bias_MLE","bias_MUE"))
##    binding  fixC ar  bias_MLE_1  bias_MLE_2  bias_MLE_3 bias_MUE_1 bias_MUE_2 bias_MUE_3
## 1:   FALSE FALSE  5 0.02337984 0.02334356 0.02434576    -0.0010    -0.0010    -0.0024
## 2:   FALSE FALSE 10 0.01441520 0.01414639 0.01574731    -0.0033    -0.0036     0.0016
## 3:   FALSE  TRUE  5 0.02337984 0.02334356 0.02434576    -0.0010    -0.0010    -0.0024
## 4:   FALSE  TRUE 10 0.01441520 0.01414639 0.01574731    -0.0033    -0.0036     0.0016
## 5:    TRUE FALSE  5 0.02242992 0.02223070 0.02338601    -0.0030    -0.0016    -0.0015
## 6:    TRUE FALSE 10 0.01296960 0.01305850 0.01413865    -0.0023    -0.0017    -0.0042
## 7:    TRUE  TRUE  5 0.02242992 0.02223070 0.02338601    -0.0030    -0.0016    -0.0015
## 8:    TRUE  TRUE 10 0.01296960 0.01305850 0.01413865    -0.0023    -0.0017    -0.0042

## ** distribution of the estimates
dt.estimate <- res2stage[decision %in% c("futility","efficacy"),]
dt.estimate[, stage := factor(stage, 1:2, c("interim","final"))]
dt.estimate[, method := paste0("method ",method)]
gg.estimate <- ggplot(dt.estimate, aes(x = estimate_MUE, fill = stage, group = stage))
## gg.estimate <- gg.estimate + geom_histogram(aes(y = ..density..), position = position_dodge()) + facet_grid(scenario~method)
gg.estimate <- gg.estimate + geom_density(alpha=0.25) + facet_grid(scenario~method)
gg.estimate

## ** special cases
ftable(reason = res2stage$reason, method = res2stage$method, scenario = res2stage$scenario)
##                               scenario    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18
## reason                 method                                                                                                   
## decreasing information 1                  0    0    1    1    0    0    1    1    0    0    0    0    0    0    0    0    0    0
##                        2                  0    0    1    1    0    0    1    1    0    0    0    0    0    0    0    0    0    0
##                        3                  0    0    1    1    0    0    1    1    0    0    0    0    0    0    0    0    0    0
## efficacy               1               3740   77 3559   67 3740   77 3559   67 3815   78 3674   67 3815   78 3674   67 3396   74
##                        2               3729   82 3554   68 3729   82 3554   68 3816   78 3677   67 3816   78 3677   67 3365   72
##                        3               4137  107 3712   83 4137  107 3712   83 4238  112 3788   83 4238  112 3788   83 3528   80
## futility               1                646 7086  603 6922  646 7086  603 6922  604 7126  571 6920  604 7126  571 6920  535 6748
##                        2                658 7120  611 6904  658 7120  611 6904  628 7180  573 6925  628 7180  573 6925  520 6742
##                        3                560 6843  579 6822  560 6843  579 6822  514 6870  535 6837  514 6870  535 6837  496 6642
## Imax reached           1                  1    1    0    0    1    1    0    0    0    0    0    0    0    0    0    0    0    0
##                        2                  1    1    0    0    1    1    0    0    0    0    0    0    0    0    0    0    0    0
##                        3                  1    1    0    0    1    1    0    0    0    0    0    0    0    0    0    0    0    0
## no boundary crossed    1               5613 2836 5838 3011 5613 2836 5838 3011 5581 2796 5755 3013 5581 2796 5755 3013 6069 3178
##                        2               5612 2797 5835 3028 5612 2797 5835 3028 5556 2742 5750 3008 5556 2742 5750 3008 6115 3186
##                        3               5302 3049 5709 3095 5302 3049 5709 3095 5248 3018 5677 3080 5248 3018 5677 3080 5976 3278

## ** reversal
res2stage.reversal <- res2stage[missing == TRUE, .(N = .N,
                                                   futility2efficacy = (stage[1] == 1)*(reason[1] == "futility")*(stage[2] == 1)*(decision[2] == "efficacy"),
                                                   efficacy2futility = (stage[1] == 1)*(reason[1] == "efficacy")*(stage[2] == 1)*(decision[2] == "futility")),
                                by = c("method","seed","binding","fixC","ar","hypo")]
res2stage.reversal[is.na(futility2efficacy), futility2efficacy := 0]
res2stage.reversal[is.na(efficacy2futility), efficacy2futility := 0]


res2stageS.reversal <- res2stage.reversal[, .(N = .N, fu2eff = 100*mean(futility2efficacy), eff2fu = 100*mean(efficacy2futility)),
                                          by = c("method","binding","fixC","ar","hypo")]
dcast(res2stageS.reversal, hypo + ar + binding + fixC ~ method, value.var = c("fu2eff","eff2fu"))
##      hypo ar binding  fixC fu2eff_1 fu2eff_2 fu2eff_3 eff2fu_1 eff2fu_2 eff2fu_3
##  1: power  5   FALSE FALSE     0.04     0.04     0.00     0.03     0.03     0.51
##  2: power  5   FALSE  TRUE     0.04     0.04     0.00     0.03     0.03     0.51
##  3: power  5    TRUE FALSE     0.06     0.08     0.02     0.05     0.07     0.65
##  4: power  5    TRUE  TRUE     0.06     0.08     0.02     0.05     0.07     0.65
##  5: power 10   FALSE FALSE     0.35     0.38     0.05     0.18     0.21     0.96
##  6: power 10   FALSE  TRUE     0.35     0.38     0.05     0.18     0.21     0.96
##  7: power 10    TRUE FALSE     0.57     0.57     0.13     0.15     0.20     1.06
##  8: power 10    TRUE  TRUE     0.57     0.57     0.13     0.15     0.20     1.06
##  9: typeI  5   FALSE FALSE     0.01     0.01     0.00     0.00     0.01     0.08
## 10: typeI  5   FALSE  TRUE     0.01     0.01     0.00     0.00     0.01     0.08
## 11: typeI  5    TRUE FALSE     0.02     0.02     0.00     0.01     0.03     0.15
## 12: typeI  5    TRUE  TRUE     0.02     0.02     0.00     0.01     0.03     0.15
## 13: typeI 10   FALSE FALSE     0.06     0.05     0.01     0.09     0.09     0.31
## 14: typeI 10   FALSE  TRUE     0.06     0.05     0.01     0.09     0.09     0.31
## 15: typeI 10    TRUE FALSE     0.11     0.11     0.03     0.09     0.08     0.36
## 16: typeI 10    TRUE  TRUE     0.11     0.11     0.03     0.09     0.08     0.36


res2stage[method == 3 & binding == FALSE & fixC == FALSE & ar == 10 & hypo == "typeI"& seed == 929803219, .(stage,type,decision,reason)]
##    stage  type decision   reason
## 1:     1 typeI     stop futility
## 2:     1 typeI efficacy     <NA>
## 3:     2 typeI     <NA>     <NA>

## ** frequency mismatch p-value / boundaries
res2stage.mismatchFU <- res2stage[decision=="futility",.(N = .N, mismatch = 100*mean(p.value_MUE<0.025)),
                                  by = c("method","binding","fixC","ar","hypo")]
dcast(res2stage.mismatchFU, hypo + ar + binding + fixC ~ method, value.var = "mismatch")
##      type ar binding  fixC         1         2          3
##  1: power  5   FALSE FALSE 0.4134367 0.4132231 0.56294780
##  2: power  5   FALSE  TRUE 0.4134367 0.4132231 0.56294780
##  3: power  5    TRUE FALSE 0.0000000 0.0000000 0.65228299
##  4: power  5    TRUE  TRUE 0.0000000 0.0000000 0.55583628
##  5: power 10   FALSE FALSE 2.4314537 2.4742268 1.34366925
##  6: power 10   FALSE  TRUE 2.4314537 2.4742268 1.34366925
##  7: power 10    TRUE FALSE 0.0000000 0.0000000 1.07416880
##  8: power 10    TRUE  TRUE 0.0000000 0.0000000 1.07416880
##  9: typeI  5   FALSE FALSE 0.0204897 0.0204876 0.06151323
## 10: typeI  5   FALSE  TRUE 0.0204897 0.0204876 0.06151323
## 11: typeI  5    TRUE FALSE 0.0000000 0.0000000 0.02049705
## 12: typeI  5    TRUE  TRUE 0.0000000 0.0000000 0.02049180
## 13: typeI 10   FALSE FALSE 0.1433692 0.1433398 0.09227930
## 14: typeI 10   FALSE  TRUE 0.1433692 0.1433398 0.09227930
## 15: typeI 10    TRUE FALSE 0.0000000 0.0000000 0.01024590
## 16: typeI 10    TRUE  TRUE 0.0000000 0.0000000 0.01024590
res2stage.mismatchEFF <- res2stage[decision=="efficacy",.(N = .N, mismatch = 100*mean(p.value_MUE>0.025)),
                                  by = c("method","binding","fixC","ar","hypo")]
all(res2stage.mismatchEFF$mismatch==0)
## [1] TRUE

## ** information
## nX1.interim is method depend. 

res2stage.nXinterim <- res2stage[missing == TRUE & method==1,.(N = .N, nX1 = unique(nX1.interim), nX2 = unique(nX2.interim), nX3 = unique(nX3.interim)),
                                 by = c("ar","seed","binding","fixC","hypo")]
## all(res2stage.nXinterim$N==3)
res2stageS.nXinterim <- res2stage.nXinterim[, .(N = .N, pc.all = 100*mean(nX3/nX1), pc.missing3 = 100*mean(nX2/nX1-nX3/nX1), pc.missing23 = 100*mean(1-nX2/nX1)),
                                            by = c("ar","hypo","fixC","binding")]
setkeyv(res2stageS.nXinterim,"ar")
res2stageS.nXinterim
##     ar  hypo  fixC binding     N   pc.all pc.missing3 pc.missing23
##  1:  5 power FALSE    TRUE 10000 79.53472    9.562374     10.90291
##  2:  5 typeI FALSE    TRUE 10000 79.53472    9.562374     10.90291
##  3:  5 power  TRUE    TRUE 10000 79.53472    9.562374     10.90291
##  4:  5 typeI  TRUE    TRUE 10000 79.53472    9.562374     10.90291
##  5:  5 power  TRUE   FALSE 10000 79.64196    9.449136     10.90890
##  6:  5 typeI  TRUE   FALSE 10000 79.64196    9.449136     10.90890
##  7:  5 power FALSE   FALSE 10000 79.64196    9.449136     10.90890
##  8:  5 typeI FALSE   FALSE 10000 79.64196    9.449136     10.90890
##  9: 10 power FALSE    TRUE 10000 71.60971   13.327969     15.06232
## 10: 10 typeI FALSE    TRUE 10000 71.60971   13.327969     15.06232
## 11: 10 power  TRUE    TRUE 10000 71.60971   13.327969     15.06232
## 12: 10 typeI  TRUE    TRUE 10000 71.60971   13.327969     15.06232
## 13: 10 power  TRUE   FALSE 10000 71.79364   13.168843     15.03752
## 14: 10 typeI  TRUE   FALSE 10000 71.79364   13.168843     15.03752
## 15: 10 power FALSE   FALSE 10000 71.79364   13.168843     15.03752
## 16: 10 typeI FALSE   FALSE 10000 71.79364   13.168843     15.03752

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
