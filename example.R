#Function to automatically load all R functions in a directory
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

#Source Paul's functions
sourceDir("R")
library(rpact) #Library for group sequential designs (and other)
library(gsDesign)

set.seed(1322)
x <- GenData(n=82*2,delta=1.3,ar=5)  #generate data with all data for in case trial completes
thear <- 10  #accrual rate (pt per month)
theDelta.t=0.7500001  #time in months to process data
thet <- x$d$t3[ceiling(nrow(x$d)/2) + ceiling(thear*theDelta.t)] #time point at which to do IA

#planned boundaries
theAlpha <- 0.025
theBeta <- 0.2
theDelta <- 1.5
theK <- 2

b1 <- CalcBoundaries(kMax=theK,  #max number of analyses (including final)
                     sided=1,  #one or two-sided
                     alpha=theAlpha,  #type I error
                     beta=theBeta,  #type II error
                     InfoR.i=c(0.5,1),  #planned or observed information rates
                     gammaA=2,  #rho parameter for alpha error spending function
                     gammaB=2,  #rho parameter for beta error spending function
                     method=1,  #use method 1 or 2 from paper H&J
                     delta=theDelta,  #effect that the study is powered for
                     InfoR.d=0.55)
plot(b1)

#planned sample size
n <- power.t.test(delta=1.5,sd=3,sig.level=0.025,power=0.8,alternative = "one.sided")$n/0.82
n*b1$InflationFactor

#decision at IA
xx <- SelectData(x$d,t=thet,Delta.t=theDelta.t*2)  #data at IA when deciding whether to continue recruitment

#time_halfpt_t3 <- x$d$t3[ceiling(nrow(x$d)/2)]
#57/10*2 #11.4 times two weeks to recruit half of the patients
#57/10*2+2 #13.4 times two weeks to recruit half of the patients and follow them up to day 28, timing of DBL
#57/10*2+2+1.5 #timing of interim analysis result is 14.9

#index11 <- which(x$d$t1<=13.4)
#xx <- x$d[index11,] #select patients with baseline data at time of DBL

#put data to NA for follow up measurements that are too late to be included in the DBL
#index12 <- which(xx$t2>13.4)
#index13 <- which(xx$t3>13.4)

#xx[index12,c("X2","X3")] <- NA
#xx[index13,"X3"] <- NA
#xx[index12,c("t2","t3")] <- NA
#xx[index13,"t3"] <- NA

b.interim <- update(b1, data = xx, k = 1, type = "interim")

IA <- analyzeData(xx) #Interim analysis results - can we add Z-value to output?
class(IA)


p <- PlotProgress(xx,Delta.t = theDelta.t)

head(xx)
nrow(xx) #92 patients have been recruited
sum(is.na(xx$t1))
sum(is.na(xx$t2))
sum(is.na(xx$t3))

dIA <- Decision(object = b1, k=1, type.k = "interim")   ## IA,Id=IA$getInformation["decision"]/b1$Imax,k=1,analysis="interim",Ik=c(IA$Info,b1$Imax)

#Decision at decision analysis
#index21 <- which(x$d$t1<=14.9)
index21 <- which(x$d$t1<=thet)

#xx_d <- x$d[x$d$id%in%as.character(c(1:100)),]
xx_d <- x$d[index21,]
#xx_d <- x$d[x$d$id%in%xx$id,]
IA_d <- AnalyzeData(xx_d) #Decision analysis results

dfinal <- Decision(IA_d,b1,Id=IA_d$Info/b1$Imax,k=1,analysis="decision",Ik=c(IA$Info,b1$Imax))

#final results
FinalPvalue(Id=IA_d$getInformation["decision"],
            Ik=IA$getInformation["interim"],
            ck=dfinal$details,
            lk=dIA$details["l"],
            uk=dIA$details["u"],
            kMax=2,
            estimate=IA_d$estimate)

FinalEstimate(Id=IA_d$getInformation["decision"],
              Ik=IA$getInformation["interim"],
              ck=dfinal$details,
              lk=dIA$details["l"],
              uk=dIA$details["u"],
              kMax=2,
              estimate=IA_d$estimate)

FinalCI(Id=IA_d$getInformation["decision"],
        Ik=IA$getInformation["interim"],
        ck=dfinal$details,
        lk=dIA$details["l"],
        uk=dIA$details["u"],
        kMax=2,
        estimate=IA_d$estimate)
