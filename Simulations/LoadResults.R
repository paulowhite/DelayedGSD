rm(list=ls())

#setwd("~/research/CI-Risk-Ratio-Diff/NewSim2021/SimuOutput/")

dirname <- "Simulations/brice-biostat-cluster/Results/simuMain/"

jobname <- "simuMainScenarioName"
alli <- 1:40

res <- NULL
for(i in alli){
  filename <- paste(dirname,jobname,'-', i, '.rda', sep="")
  if(file.exists(filename)) {
    load(file=filename)
    res <- rbind(res,RES)
  }
}
dim(res)
colnames(res)
summary(res)

table(res$method)
head(res)

res$final.efficacy <- res$statistic >= res$uk | res$statistic >= res$ck
res$final.efficacy[res$type%in%"interim"] <- NA

res$final.futility <- res$statistic <= res$lk | res$statistic < res$ck
res$final.futility[res$type%in%"interim"] <- NA

nsim <- length(unique(res$seed))
true_eff <- 0.8
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

discrep <- res[res$p.value_MUE>0.05 & !is.na(res$p.value_MUE),]
discrep <- discrep[!is.na(discrep$final.efficacy),]
dim(discrep)

res[res$seed%in%44942,]
#weird that decision indicates futility. If I rerun the simulation iteration myself it concludes efficacy, which is correct
#weird scenario with info at final equal to info at interim... This is an error in the code to make the output of the simulations. However, this doesn't solve the p-value issue

res[res$seed%in%79729,]
#same issue with info at final equal to info at interim...

res[res$seed%in%32877,]
#same issue with info at final equal to info at interim...
#p-value isn't even close to 0.05


#H0 rejected based on boundaries, final p-value, final MLE, final MUE, CI, was there a flip, weird case

#missing: weird cases identifier and flip identifier

#power per stage and total
#mean bias
#median bias
#CI coverage
#distribution of MLE estimate
#flipping probabilities
#how often do weird special cases occur
#is p-value always < 0.05 when H0 rejected?
#distribution of p-value