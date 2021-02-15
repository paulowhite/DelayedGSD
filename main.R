### main.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 24 2020 (14:37) 
## Version: 
## Last-Updated: Feb 12 2021 (21:25) 
##           By: Paul Blanche
##     Update #: 336
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
## res <- GenData(n=100,N.fw=4,ar=10)

## MyMissProb <- matrix(c(4,1,6,93),ncol=2,nrow=2,byrow=TRUE)
MyMissProb <- matrix(c(5,1,6,92),ncol=2,nrow=2,byrow=TRUE) # to additionnally remove 1 more because some FASFL=N
colnames(MyMissProb) <- c("V2 missing","V2 not missing")
rownames(MyMissProb) <- c("V1 missing","V1 not missing")
MyMissProb <- MyMissProb/sum(MyMissProb)

res <- GenData(n=104, 
               N.fw=2,
               rand.block=c(1,1,0,0),
               allsd=c(2.5,2.1,2.4),
               mean0=c(10,0,0),
               delta=c(0,-0.6,-0.8),
               ar=(0.86*2)*2,
               cor.01.1=-0.15,
               cor.ij.1=0.68,
               cor.0j.1=-0.27,
               seed=24082020,
               MissProb=MyMissProb)
d <-res$d
head(d)
tail(d)
dim(d)
summary(d)

# {{{ names of variable and format as in Corine's example
## Make data long format
long <- reshape(d,
                 direction='long',
                 varying=c('X2','X3'),
                 idvar='id',
                 v.names='Change')
long <- long[order(long$id),]
head(d)
head(long)
#----
dd <- long
names(dd) <- c("USUBJID","TRT01P","BASE","RANDDT","DT14","DT28","AVISIT","CHG")
dd$TRT01P <- factor(dd$TRT01P,levels=c(0,1),labels=c("Control","Active"))
dd$AVISIT <- factor(dd$AVISIT,levels=c(1,2),labels=c("Day 14","Day 28"))
# time unit in days
dd$RANDDT <- round(dd$RANDDT*14)
dd$DT14 <- round(dd$DT14*14)
dd$DT28 <- round(dd$DT28*14)
# round with 2 digits
dd$BASE <- round(dd$BASE,2)
dd$CHG <- round(dd$CHG,2)
head(dd)
# Add FASFL
IDFASFLN <- d[is.na(res$d$X2) & is.na(res$d$X3),"id"]
dd$FASFL <- factor(dd$USUBJID %in% IDFASFLN,levels=c(TRUE,FALSE),labels=c("N","Y"))
summary(dd)

# {{{ save data example
write.table(dd,file="~/research/SeqDesignDelayed/DelayedGSD/ExampleData.csv",sep=",",row.names=FALSE)
# }}}


#---------------
# similar analysis as in Corine's example
## library(nlme)
m <- gls(CHG~BASE*AVISIT + TRT01P*AVISIT,
         ## data=dd,
         data=dd[dd$FASFL=="Y",],
         correlation=corSymm(form=~1|USUBJID),
         varIdent(form=~1|AVISIT),
         method="REML",
         na.action = na.exclude)
summary(m)
# }}}


# {{{ check data generation seems fine
cat("\n",paste0("Accrual rate, chosen= ",res$input$ar,", obs. =",round(1/mean(diff(d[,"t1"])),4)," per unit of time between 2 follow-up visits."),"\n")
cat("\n",paste0("Mean in placebo, chosen= ",paste(res$input$mean0,collapse="-"),", obs. =",paste(round(colMeans(res$d[res$d$Z==0,paste0("X",(1:(res$input$N.fw+1)))],na.rm=TRUE),4),collapse=";")," per unit of time between 2 follow-up visits."),"\n")
cat("\n",paste0("Difference in Means, chosen= ",paste(res$input$delta,collapse="-"),", obs. =",
       paste(round(colMeans(res$d[res$d$Z==1,paste0("X",(1:(res$input$N.fw+1)))],na.rm=TRUE)-colMeans(res$d[res$d$Z==0,paste0("X",(1:(res$input$N.fw+1)))],na.rm=TRUE),4),collapse=";")
      ," ."),"\n")
cat("\n",paste0("All sd, chosen= ",paste(res$input$allsd,collapse="-"),", obs. =",paste(round(apply(res$d[res$d$Z==0,paste0("X",(1:(res$input$N.fw+1)))],2,sd,na.rm=TRUE),4),collapse="-"),"."),"\n")
cat("\n",paste0("Corr(X1,X2), chosen= ",res$input$cor.01.1,", obs. =",round(cor(res$d$X1,res$d$X2,use="pairwise.complete.obs"),4)),"\n")
cat("\n",paste0("Corr(X2,X3), chosen= ",res$input$cor.ij.1,", obs. =",round(cor(res$d$X2,res$d$X3,use="pairwise.complete.obs"),4)),"\n")
cat("\n",paste0("Corr(X1,X3),  obs. =",round(cor(res$d$X1,res$d$X3,use="pairwise.complete.obs"),4)),"\n")


cat("\n",paste0("Missing at V1 and V2, chosen= ",round(res$input$MissProb[1,1]*100,2),"% , obs. =",round(mean(is.na(res$d$X2) & is.na(res$d$X3))*100,2)),"% \n")
cat("\n",paste0("Missing at V1 and not V2, chosen= ",round(res$input$MissProb[1,2]*100,2),"% , obs. =",round(mean(is.na(res$d$X2) & !is.na(res$d$X3))*100,2)),"% \n")
cat("\n",paste0("Missing at V2 and not V1, chosen= ",round(res$input$MissProb[2,1]*100,2),"% , obs. =",round(mean(!is.na(res$d$X2) & is.na(res$d$X3))*100,2)),"% \n")
# }}}


d[is.na(res$d$X2) & is.na(res$d$X3),]



#
#


#----- select when half of the subjects have one follow-up measuement---
# which depends on accrual rate (ar) and time to process data (Delta.t)
theDelta.t <- 0.500001
thear <- 10
thet <- d$t2[ceiling(nrow(d)/2) + ceiling(thear*theDelta.t)]

# plot progress in including patients and collected data
PlotProgress(d,at=thet,Delta.t=0.5)

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

coefzoom <- 2
## png(filename = "~/research/SeqDesignDelayed/BoundaryComparison.png",
    ## width = 480*coefzoom, height = 480*coefzoom, units = "px")
par(mfrow=c(2,2))
PlotBoundaries(PlannedB1F,type="Z",Itype="rate")
title("Method 1, c value at decision allowed too low")
PlotBoundaries(PlannedB2F,type="Z",Itype="rate")
title("Method 2, c value at decision allowed too low")
PlotBoundaries(PlannedB1T,type="Z",Itype="rate")
title("Method 1, c value at decision NOT allowed too low")
PlotBoundaries(PlannedB2T,type="Z",Itype="rate")
title("Method 2, c value at decision NOT allowed too low")
## dev.off()

#----------------------------------------------------------------------
### main.R ends here

