rm(list=ls())
res <- readRDS("Simulations/brice-biostat-cluster/Results-built/res2stage.rds")
dim(res)

res$final.efficacy <- res$statistic >= res$uk | res$statistic >= res$ck
res$final.efficacy[res$type%in%"interim"] <- NA

res$final.futility <- res$statistic <= res$lk | res$statistic < res$ck
res$final.futility[res$type%in%"interim"] <- NA

result <- NULL

for(scen in 1:length(unique(res$scenario))){
  print(scen)
  res_temp2 <- res[res$scenario%in%scen,]
  nsim <- length(unique(res_temp2$seed))
  
  true_eff <- ifelse(unique(res_temp2$hypo)%in%"power",0.6,0)
  
  for(m in c(1:3)){
    print(m)
    res_temp <- res_temp2[res_temp2$method%in%m,]
    result <- rbind(result,c("scenario"=scen,
                             "missing"=unique(res_temp$missing),
                             "binding"=unique(res_temp$binding),
                             "fixC"=unique(res_temp$fixC),
                             "ar"=unique(res_temp$ar),
                             "hypo"=unique(res_temp$hypo),
                           "method"=m,
                           "power_bnds"=round(sum(res_temp$final.efficacy,na.rm=T)/nsim,4),
                           "power_pval"=round(sum(res_temp$p.value_MUE<0.025,na.rm=T)/nsim,4),
                           "discrep_pval_bnds"=round((sum(res_temp$p.value_MUE<0.025 & is.na(res_temp$final.efficacy),na.rm=T)+sum(res_temp$p.value_MUE>0.025 & res_temp$final.efficacy,na.rm=T))/nsim,4),
                           "bias_MLE"=round(sum(res_temp$estimate_ML[res_temp$type%in%c("decision","final")],na.rm=T)/nsim-true_eff,4),
                           "bias_MUE"=round(sum(res_temp$estimate_MUE[res_temp$type%in%c("decision","final")]>true_eff,na.rm=T)/nsim,4),
                           "CI_coverage"=round(sum(true_eff>=res_temp$lower_MUE[res_temp$type%in%c("decision","final")] & true_eff<=res_temp$upper_MUE[res_temp$type%in%c("decision","final")],na.rm=T)/nsim,4)
    ))
    #hist(res_temp$p.value_MUE)
  }
}

save(result,file="Simulations/ResTwoStage24102022.Rdata")

discrep <- res[res$p.value_MUE>0.025 & !is.na(res$p.value_MUE),]
discrep <- discrep[!is.na(discrep$final.efficacy),]
dim(discrep)

discrep <- res[res$p.value_MUE<0.025 & !is.na(res$p.value_MUE),]
discrep <- discrep[is.na(discrep$final.efficacy),]
dim(discrep)