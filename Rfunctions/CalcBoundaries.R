CalcBoundaries <- function(kMax=2,  #max number of analyses (including final)
                           sided=1,  #one or two-sided
                           alpha=0.025,  #type I error
                           beta=0.2,  #type II error
                           informationRates=c(0.5,1),  #planned or observed information rates
                           gammaA=2,  #rho parameter for alpha error spending function
                           gammaB=2,  #rho parameter for beta error spending function
                           method=1,  #use method 1 or 2 from paper H&J
                           cNotBelowFixedc=FALSE, # whether the value c at the decision analysis can be below that of a fixed sample test (H & J page 10)
                           delta=1.5, #effect that the study is powered for
                           Id=0.55,   #(expected) information ratio at each decision analysis
                           bindingFutility=TRUE,  #whether the futility stopping rule is binding
                           Trace=TRUE){  # whether to print some messages
                            
  require(rpact)
  
  if(sided!=1){
    stop("Function cannot handle two-sided tests yet")
  }
  
  if(sum(Id < informationRates[-length(informationRates)])>0){
    stop("Information at decision analysis should not be smaller than information at interim analysis")
  }
  
  
  if(Trace){message("In CalcBoundaries, the method assumes that positive effects are good")}
  
  #it looks like Rpact cannot handle binding futility rules... I have switched to using gsDesign
  #StandardDesign <- getDesignGroupSequential(kMax=kMax, sided = sided, alpha = alpha, beta = beta,
  #                                          informationRates = informationRates,
  #                                          typeOfDesign="asKD",
  #                                          typeBetaSpending="bsKD",gammaA=gammaA,gammaB=gammaB, 
  #                                          bindingFutility = bindingFutility)
  
  if(bindingFutility){
    StandardDesign <- gsDesign(k=kMax, test.type=3,alpha=alpha,beta=beta,
                               timing=informationRates,
                               sfu=sfPower,sfupar=gammaA,
                               sfl=sfPower,sflpar=gammaB)
  } else {
    StandardDesign <- gsDesign(k=kMax, test.type=4,alpha=alpha,beta=beta,
                               timing=informationRates,
                               sfu=sfPower,sfupar=gammaA,
                               sfl=sfPower,sflpar=gammaB)
  }
  
  
  #R <- getDesignCharacteristics(StandardDesign)$inflationFactor
  R <- StandardDesign$n.I[kMax]
  Imax <- ((qnorm(1-alpha)+qnorm(1-beta))/delta)^2*R
  Id <- Id*Imax
  
  if(max(Id) >= Imax){
    stop("Function cannot handle Id >= Imax yet")
  }
    
  #uk <- StandardDesign$criticalValues
  #lk <- c(StandardDesign$futilityBounds,uk[kMax])
  uk <- StandardDesign$upper$bound
  lk <- StandardDesign$lower$bound
  ck <- rep(0,kMax-1)
  
  #browser()
  
  if(method==1){
    
      for(i in 1:(kMax-1)){
          ck[i] <- Method1(uk=uk[1:i],
                           lk=lk[1:i],
                           Ik=(informationRates*Imax)[1:i],
                           Id=Id[i],
                           Imax=Imax,
                           cMin=ifelse(cNotBelowFixedc,qnorm(1-alpha),-Inf),
                           bindingFutility=bindingFutility)
      }
    
  } else if(method==2){
    
      for(i in 1:(kMax-1)){
          ###THIS GIVES ERRORS BECAUSE IT HAS TROUBLE WITH THE BOUNDARIES FOR THE UNIROOT FUNCTIONS IN METHOD1 AND 2
          delayedBnds <- Method2(uk=uk[1:i],
                                 lk=lk[1:i],
                                 Ik=(informationRates*Imax)[1:i],
                                 Id=Id[i],
                                 Imax=Imax,
                                 delta=delta,
                                 cMin=ifelse(cNotBelowFixedc,qnorm(1-alpha),-Inf),
                                 bindingFutility=bindingFutility)
          ck[i] <- delayedBnds$critval
          lk[i] <- delayedBnds$lkd
      }
    
  } else {
    stop(("Please specify method=1 or method=2"))
  }
  
  list("uk"=uk,"lk"=lk,"ck"=ck,"Ik"=informationRates*Imax,"Id"=Id,"Imax"=Imax,
       "alpha"=alpha,kMax=kMax,sided=sided,beta=beta,gammaA=gammaA,gammaB=gammaB,method=method,delta=delta,cNotBelowFixedc=cNotBelowFixedc,bindingFutility=bindingFutility)
}

