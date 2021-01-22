#Function to automatically load all R functions in a directory
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

#Source Paul's functions
sourceDir("Rfunctions")
library(rpact) #Library for group sequential designs (and other)

x <- GenData(n=68*2)  #generate data with all data for in case trial completes
thear <- 10  #accrual rate (pt per month)
theDelta.t=0.500001  #time in months to process data
thet <- x$d$t2[ceiling(nrow(x$d)/2) + ceiling(thear*theDelta.t)] #time point at which to do IA
xx <- SelectData(x$d,t=thet,Delta.t=theDelta.t)  #data at IA when deciding whether to continue recruitment
IA <- AnalyzeData(xx) #Interim analysis results

#planned boundaries
theAlpha <- 0.025
theBeta <- 0.2
theDelta <- 1.5
theK <- 2

b1 <- CalcBoundaries(kMax=theK,  #max number of analyses (including final)
                     sided=1,  #one or two-sided
                     alpha=theAlpha,  #type I error
                     beta=theBeta,  #type II error
                     informationRates=c(0.5,1),  #planned or observed information rates
                     gammaA=2,  #rho parameter for alpha error spending function
                     gammaB=2,  #rho parameter for beta error spending function
                     method=1,  #use method 1 or 2 from paper H&J
                     delta=theDelta,  #effect that the study is powered for
                     Id=0.55)

b2 <- CalcBoundaries(kMax=theK,  #max number of analyses (including final)
                     sided=1,  #one or two-sided
                     alpha=theAlpha,  #type I error
                     beta=theBeta,  #type II error
                     informationRates=c(0.5,1),  #planned or observed information rates
                     gammaA=2,  #rho parameter for alpha error spending function
                     gammaB=2,  #rho parameter for beta error spending function
                     method=2,  #use method 1 or 2 from paper H&J
                     delta=theDelta,  #effect that the study is powered for
                     Id=0.55)

#boundary that used to produce an error message (not anymore thanks to Brice!). It still produces some warnings though...
b22 <- CalcBoundaries(kMax=3,  #max number of analyses (including final)
                     sided=1,  #one or two-sided
                     alpha=theAlpha,  #type I error
                     beta=theBeta,  #type II error
                     informationRates=c(0.6,0.8,1),  #planned or observed information rates
                     gammaA=2,  #rho parameter for alpha error spending function
                     gammaB=2,  #rho parameter for beta error spending function
                     method=2,  #use method 1 or 2 from paper H&J
                     delta=theDelta,  #effect that the study is powered for
                     Id=c(0.61,0.81))

PlotBoundaries(b2)

#plot planned boundaries
par(mfrow=c(2,3))
  PlotBoundaries(b1,type="Z",main="Method 1")
  PlotBoundaries(b1,type="P",main="Method 1")
  PlotBoundaries(b1,type="E",main="Method 1")
  PlotBoundaries(b2,type="Z",main="Method 2")
  PlotBoundaries(b2,type="P",main="Method 2")
  PlotBoundaries(b2,type="E",main="Method 2")


#Decision regarding recruitment at interim analysis
Decision(IA,b1,Id=0.55,k=1,analysis="interim",Ik=c(IA$Info,b1$Imax))  
Decision(IA,b2,Id=0.55,k=1,analysis="interim",Ik=c(IA$Info,b2$Imax))    

#Decision at decision analysis
xx_d <- x$d[x$d$id%in%xx$id,]
IA_d <- AnalyzeData(xx_d) #Decision analysis results

Decision(IA_d,b1,Id=IA_d$Info/b1$Imax,k=1,analysis="decision",Ik=c(IA$Info,b1$Imax))
Decision(IA_d,b2,Id=IA_d$Info/b2$Imax,k=1,analysis="decision",Ik=c(IA$Info,b2$Imax))  

#Decision at final analysis
FA <- AnalyzeData(x$d)
Decision(FA,b1,Id=IA_d$Info/b1$Imax,k=2,analysis="final",Ik=c(IA$Info,FA$Info))  
Decision(FA,b2,Id=IA_d$Info/b2$Imax,k=2,analysis="final",Ik=c(IA$Info,FA$Info))    


######################### OLD ####################################################
plannedDesign <- getDesignGroupSequential(kMax=2, sided = 1, alpha = theAlpha, beta = theBeta,
                                          informationRates = c(0.5,1),
                                          typeOfDesign="asKD",
                                          typeBetaSpending="bsKD",gammaA=2,gammaB=2)

#function to update a planned GSD with observed information rates (ignoring delayed outcomes)
updateDesign <- function(plannedDesign,informationRates){
  getDesignGroupSequential(kMax=plannedDesign$kMax, sided = plannedDesign$sided, alpha = plannedDesign$alpha, beta = plannedDesign$beta,
                           informationRates = informationRates,
                           typeOfDesign=plannedDesign$typeOfDesign,
                           typeBetaSpending=plannedDesign$typeBetaSpending,gammaA=plannedDesign$gammaA,gammaB=plannedDesign$gammaB)
}

#calculate critval c using method 1:
method1 <- function(uk,  #upper bounds for all analyses up to and including current stage k
                    lk,  #lower bounds for all analyses up to and including current stage k
                    Ik,  #Information for all analyses up to and including current stage k
                    Id){  #Observed information at decision analysis k
  k <- length(uk)
  Ik <- c(Ik,Id)
  sigmaZk <- diag(1,k+1)
  for(i in 1:(k+1)){
    for(j in i:(k+1)){
      sigmaZk[i,j] <- sqrt(Ik[i]/Ik[j])
      sigmaZk[j,i] <- sqrt(Ik[i]/Ik[j])
    }
  }
  
  c <- uniroot(function(x){pmvnorm(lower = c(lk[0:(k-1)],uk[k],-Inf),
                                   upper = c(uk[0:(k-1)],Inf,x),
                                   mean=rep(0,k+1),
                                   sigma= sigmaZk) - 
      pmvnorm(lower = c(lk[0:(k-1)],-Inf,x),
              upper = c(uk[0:(k-1)],lk[k],Inf),
              mean=rep(0,k+1),
              sigma= sigmaZk)},
      lower = lk[k],
      upper = uk[k])$root
  c
}

# rho-family spending functions (Kim-DeMets) for alpha and beta 
g <- function(I,rho,beta,Imax){
  beta*min(1,(I/Imax)^rho)
}

method2 <- function(uk,  #upper bounds for all analyses up to and including current stage k
                    lk,  #lower bounds for all analyses up to and including current stage k
                    Ik,  #Information for all analyses up to and including current stage k
                    Id,  #Expected information at decision analysis
                    Imax, #planned maximum information
                    delta){ #planned effect
  
  k <- length(uk)
  Ik <- c(Ik,Id)
  sigmaZk <- diag(1,k+1)
  for(i in 1:(k+1)){
    for(j in i:(k+1)){
      sigmaZk[i,j] <- sqrt(Ik[i]/Ik[j])
      sigmaZk[j,i] <- sqrt(Ik[i]/Ik[j])
    }
  }
  
  theta <- delta*sqrt(Ik)
  
  sol <- uniroot(function(x){
    
    c <- method1(lk=c(lk[0:(k-1)],x),uk=uk,Ik=Ik[1:k],Id=Id)
    
    pmvnorm(lower = c(lk[0:(k-1)],uk[k],-Inf),
            upper = c(uk[0:(k-1)],Inf,c),
            mean=theta,
            sigma= sigmaZk) + 
      pmvnorm(lower = c(lk[0:(k-1)],-Inf,-Inf),
              upper = c(uk[0:(k-1)],x,c),
              mean=theta,
              sigma= sigmaZk) - (g(I=Id,Imax=Imax,beta=0.2,rho=2)-g(I=Ik[k],Imax=Imax,beta=0.2,rho=2))},
    lower = -10, #needs something better
    upper = uk[k]*0.999)$root #can't be exactly uk[k] since upper bound for method1 function and then the lower and upper bound would be equal
  
  c <- method1(lk=c(lk[0:(k-1)],sol),uk=uk,Ik=Ik[1:k],Id=Id)
  
  list(lkd=sol,critval=c)
}

IA_decision <- function(resultsIA,  #results from IA number k
                     k=1,  #which IA
                     plannedDesign,  #the planned boundaries
                     informationRates,  #vector of length K including the information rates as observed up to interim analsis k and as expected from analysis k+1 to K
                     plannedImax,  #planned maximum information
                     method=1,  #which method should be used to calculated the boundaries
                     Id=NULL,  #observed information at the decision analysis (method 1), or expected information at the decision analysis (method 2)
                     analysis="interim",  #is it an interim or a decision analysis
                     delta=NULL){  #the effect the trial is powered for
  currentDesign <- updateDesign(plannedDesign,informationRates)
  
  if(method==1){
    if(analysis=="interim"){
      if(resultsIA$estimate/resultsIA$se > currentDesign$criticalValues[k]){
          "stop_rec_eff"
        } else if(resultsIA$estimate/resultsIA$se < currentDesign$futilityBounds[k]){
          "stop_rec_fut"
        } else {
          "continue"
        }
    } else if(analysis=="decision"){
      c <- method1(uk=currentDesign$criticalValues[1:k],lk=currentDesign$futilityBounds[1:k],Ik=informationRates*plannedImax[1:k],Id=resultsIA$Info)
      
      if(resultsIA$estimate/resultsIA$se > c){
        "efficacy"
      } else {
        "futility"
      }
    } else {
      warning("please specify analysis=interim or analysis=decision")
    }
  } else if(method==2){
    delayedBnds <- method2(uk=currentDesign$criticalValues[1:k],
                           lk=currentDesign$futilityBounds[1:k],
                           Ik=informationRates*plannedImax[1:k],
                           Id=Id,
                           Imax=plannedImax,
                           delta=delta)
    
    if(analysis=="interim"){
      if(resultsIA$estimate/resultsIA$se > currentDesign$criticalValues[k]){
        "stop_rec_eff"
      } else if(resultsIA$estimate/resultsIA$se < delayedBnds$lkd){
        "stop_rec_fut"
      } else {
        "continue"
      }
    } else if (analysis=="decision"){
      if(resultsIA$estimate/resultsIA$se > delayedBnds$critval){
        "stop_rec_eff"
      } else {
        "futility"
        }
    } else {
      warning("please specify analysis=interim or analysis=decision")
    }
    
  } else {
    warning("Please specify method=1 or method=2")
  }
}

dchar <- getDesignCharacteristics(plannedDesign) #Design characteristics, amongst others the inflation factor
thePlannedImax <- ((qnorm(1-theAlpha)+qnorm(1-theBeta))/theDelta)^2 * dchar$inflationFactor  #calculate required max info to obtain power for planned effect
informationRates <- c(IA$Info/Imax,1) #information rates (actual for IA and planned for final)

#results at interim analysis using method 1:
IA_decision(resultsIA=IA,
            k=1,
            plannedDesign=plannedDesign,
            informationRates=c(IA$Info,thePlannedImax)/thePlannedImax,
            plannedImax=thePlannedImax)

#results at decision analysis:


#need a final analysis function using the updated design as in the rpact tutorial


##############OLD plotting
#trying to plot stuf. Might be nice to create a plot function that can plot the boundaries
Irates <- c(0.5,1)
Idrates <- 0.55
plot(Irates,b2$uk,type="l",ylim=c(-1,3),col="blue",xlab="Information fraction",ylab="Stopping boundary (Z-statistic)")
lines(Irates,b2$lk,col="blue")
points(Idrates,b2$ck,col="blue")
points(Irates,b2$lk,col="blue")
points(Irates,b2$uk,col="blue")
lines(Irates,b1$uk,col="red")
lines(Irates,b1$lk,col="red")
points(Irates,b1$uk,col="red")
points(Irates,b1$lk,col="red")
points(Idrates,b1$ck,col="red")
legend("bottomright",c("Method 1","Method 2"),lty=c(1,1),col=c("red","blue"),pch=c(1,1))


f <- plot(StandardDesign)
f+geom_point(aes(x=(Id/Imax)[1], y=ck[1]), colour="purple",size=3)+
  geom_point(aes(x=(Id/Imax)[2], y=ck[2]), colour="purple",size=3)+
  geom_point(aes(x=(Id/Imax)[1], y=ck2[1]), colour="green",size=3)+
  geom_point(aes(x=(Id/Imax)[2], y=ck2[2]), colour="green",size=3)+
  geom_point(aes(x=(Ik/Imax)[1], y=lk[1]), colour="green",size=3)+
  geom_point(aes(x=(Ik/Imax)[2], y=lk[2]), colour="green",size=3)

geom_point(aes(x=Ik/Imax, y=c2$lkd), colour="lightblue",size=3)