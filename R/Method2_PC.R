#' @title Calculate boundaries for a group sequential design using Method 2
#' @description Calculate boundaries for a group sequential design with delayed endpoints based on planned and/or observed information using an error spending approach with a binding futility boundary.
#' 
#' 
#' 
#' @param rho_alpha rho parameter of the rho-family spending functions (Kim-DeMets) for alpha
#' @param rho_beta rho parameter of the rho-family spending functions (Kim-DeMets) for beta
#' @param alpha type I error
#' @param beta type II error
#' @param kMax max number of analyses (including final)
#' @param Info.max maximum information needed for given beta (type II-error), delta (expected difference), alpha (type I-error) and Kmax. It can be given if it is known. Otherwise it is computed from the  values given for alpha, beta, delta and Kmax.
#' @param InfoR.i Expected or observed (wherever possible) information rates at the interim and final analyses 1:Kmax
#' @param InfoR.d Expected or observed information rates at all potential decision analyses 1:(Kmax-1)
#' @param delta expected effect under the alternative (should be on the scale of the test statistc for which If and Info.max relate to one over the variance, e.g. delta=expected log(Hazard ratio))
#' @param abseps tolerance for precision when finding roots or computing integrals
#' @param direction greater is for Ho= theta > 0, "smaller" is for Ho= theta < 0 (note that in Jennison and Turnbull's book chapter (2013) they consider smaller)
#' @param Trace Used only if Info.max=NULL. Whether to print informations to follow the progression of the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax).
#' @param nWhileMax Used only if Info.max=NULL. Maximum number of steps in the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax)
#' @param toldiff Used only if Info.max=NULL. Maximum tolerated difference between lower and upper bounds at anaylis Kmax (which souhld be zero), in the root finding algorithm, to find the value of Info.max
#' @param tolcoef Used only if Info.max=NULL. Maximum tolerated difference before stopping the search (in the root finding algorithm), between two successive values for the multiplier coeficient 'coef' such that Info.max=coef*If  (some values for coef are given in Table 7.6 page 164 Jennison's book. The value of "If" (which stands for Information for Fixed design) corresponds to the information we need if Kmax=1)
#' @param mycoefMax Used only if Info.max=NULL. Upper limit of the interval of values in which we search for the multiplier coeficient 'coef' such that Info.max=coef*If (in the root finding algorithm).
#' @param mycoefL Used only if Info.max=NULL. Lower limit of the interval (see mycoefMax)
#' @param myseed seed for producing reproducible results. Because we call functions which are based on Monte-Carlo compuation (pmvnorm)
#' @param cMin minimun possible value c for the decision analysis, typically that for a fixed sample test (H & J page 10)
#' @param sided one or two sided
#'  
#'
#' @examples
#'
#' Example to check that code matches
#' b1 <- CalcBoundaries(kMax=2,  #max number of analyses (including final)
#'                     sided=1,  #one or two-sided
#'                     alpha=0.025,  #type I error
#'                     beta=0.2,  #type II error
#'                     InfoR.i=c(0.6,1),  #planned information rates
#'                     rho_alpha=2,  #rho parameter for alpha error spending function
#'                     rho_beta=2,  #rho parameter for beta error spending function
#'                     method=1,  #use method 1 or 2 from paper H&J
#'                     delta=1.5,  #effect that the study is powered for
#'                     InfoR.d=0.65,
#'                     bindingFutility=TRUE)
#'
#' b12 <- Method2_PC(Kmax=2,Info.max=b1$Info.max,delta=1.5,alpha=0.025,InfoR.i=c(0.6,1),InfoR.d=0.65)
#' 
#' 
#' b1FT <- CalcBoundaries(kMax=2,  #max number of analyses (including final)
#'                     sided=1,  #one or two-sided
#'                     alpha=0.025,  #type I error
#'                     beta=0.2,  #type II error
#'                     InfoR.i=c(0.6,1),  #planned information rates
#'                     rho_alpha=2,  #rho parameter for alpha error spending function
#'                     rho_beta=2,  #rho parameter for beta error spending function
#'                     method=1,  #use method 1 or 2 from paper H&J
#'                     delta=1.5,  #effect that the study is powered for
#'                     InfoR.d=0.65,
#'                     bindingFutility=TRUE)
#'
#' b12FT <- Method2_PC(Kmax=2,Info.max=b1FT$Info.max,delta=1.5,binding=TRUE,alpha=0.025,InfoR.i=c(0.6,1),InfoR.d=0.65)
#' 
#' b12FTNoImax <- Method2_PC(Kmax=2,Info.max=NULL,delta=1.5,binding=TRUE,alpha=0.025,InfoR.i=c(0.6,1),InfoR.d=0.65)
#' 
#' 
#' 
#' all.equal(b1$uk, b12$boundaries[,"b.k"])
#' all.equal(b1$lk, b12$boundaries[,"a.k"])
#' 
#' b13 <- Method2_PC(Kmax=2,delta=1.5,alpha=0.025,Trace=T,InfoR.i=c(0.6,1))
#' 
#' 
#' to reproduce bounds from CJ DSBS course slide 106 REQUIRES NON-BINDING FUTILITY
#'               
#' bCJ <- Method2_PC(rho_alpha=1.345,
#'            rho_beta=1.345,
#'            alpha=0.025,
#'            beta=0.1,
#'            Kmax=3,
#'            Info.max=12,
#'            InfoR.i=c(3.5,6.75,12)/12,
#'            InfoR.d=c(5.5,8.75)/12,
#'            delta=1,  
#'            abseps = 1e-06, 
#'            direction="smaller",
#'            sided=1,
#'            cMin=1.96
#'            )




## * Method2_PC (code)
Method2_PC <- function(rho_alpha=2,          # rho parameter of the rho-family spending functions (Kim-DeMets) for alpha
                         rho_beta=2,           # rho parameter of the rho-family spending functions (Kim-DeMets) for beta
                         alpha=0.025,          # Type-I error (overall)
                         beta=0.2,             # Type-II error (overall)
                         Kmax,                 # number of planned analyses (including the final analysis)
                         Info.max=NULL,        # Info.max, i.e. maximum information needed for given beta (type II-error), delta (expected difference), alpha (type I-error) and Kmax. It can be given if it is known. Otherwise it is computed from the  values given for alpha, beta, delta and Kmax.
                         InfoR.i=NULL,         # Expected or observed (wherever possible) information rates at the interim and final analyses 1:Kmax
                         InfoR.d=NULL,         # Expected or observed information rates at all potential decision analyses 1:(Kmax-1)
                         delta=0,              # expected effect under the alternative (should be on the scale of the test statistc for which If and Info.max relate to one over the variance, e.g. delta=expected log(Hazard ratio))
                         abseps = 1e-06,       # tolerance for precision when finding roots or computing integrals
                         direction="smaller",  # greater is for Ho= theta > 0, "smaller" is for Ho= theta < 0 (note that in Jennison and Turnbull's book chapter (2013) they consider smaller)
                         Trace=FALSE,          # Used only if Info.max=NULL. Whether to print informations to follow the progression of the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax).
                         nWhileMax=30,         # Used only if Info.max=NULL. Maximum number of steps in the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax)
                         toldiff= 1e-05,       # Used only if Info.max=NULL. Maximum tolerated difference between lower and upper bounds at anaylis Kmax (which souhld be zero), in the root finding algorithm, to find the value of Info.max
                         tolcoef= 1e-04,       # Used only if Info.max=NULL. Maximum tolerated difference before stopping the search (in the root finding algorithm), between two successive values for the multiplier coeficient 'coef' such that Info.max=coef*If  (some values for coef are given in Table 7.6 page 164 Jennison's book. The value of "If" (which stands for Information for Fixed design) corresponds to the information we need if Kmax=1)
                         mycoefMax= 1.2,       # Used only if Info.max=NULL. Upper limit of the interval of values in which we search for the multiplier coeficient 'coef' such that Info.max=coef*If (in the root finding algorithm).
                         mycoefL=1,            # Used only if Info.max=NULL. Lower limit of the interval (see mycoefMax)
                         myseed=2902,           # seed for producing reproducible results. Because we call functions which are based on Monte-Carlo compuation (pmvnorm)
                         sided=1,              # one or two sided
                         cMin=-Inf           # minimun possible value c for the decision analysis, typically that for a fixed sample test (H & J page 10)
){
  ## {{{ set seed
  old <- .Random.seed # to save the current seed
  set.seed(myseed)
  ## }}}
  ## {{{ preliminaries
  mycoef <- NULL # initialize as needed for output
  lk <- rep(-Inf,Kmax)
  uk <- rep(Inf,Kmax) 
  ck <- rep(NA,Kmax)
  
  if( (direction=="smaller" & delta<0) | (direction=="greater" & delta>0)){
    stop("The values given for arguments direction and delta are inconsistent.\n When direction=smaller, delta should be positive.\n When direction=greater, delta should be negative.")
  }
  
  if(direction=="greater"){
    delta <- -delta
  }else if(!direction%in%c("greater","smaller")){
    stop("direction should be either greater or smaller")
  }
  
  # initialize
  thealpha <- rep(0,Kmax)  # alpha spent up to step k
  thebeta <- rep(0,Kmax)   # beta spent up to step k
  IncAlpha <- rep(0,Kmax)  # alpha spent at step k
  IncBeta <- rep(0,Kmax)   # beta spent at step k     
  # compute variance-covariance matrix of vector (Z_1,...,Z_k)
  sigmaZk <- diag(1,Kmax)
  for(i in 1:Kmax){
    for(j in i:Kmax){
      sigmaZk[i,j] <- sqrt(InfoR.i[i]/InfoR.i[j])
      sigmaZk[j,i] <- sqrt(InfoR.i[i]/InfoR.i[j])
    }
  }
  # compute If (see Jennison book page 87
  If <- (qnorm(1-alpha)+qnorm(1-beta))^2/delta^2
  if(Trace){
    cat("\n If computed as =",If,"\n")
  }
  ## }}}
  if(is.null(Info.max)){
    #browser()
    ## {{{ Compute Info.max from If and the other arguments (Recursive calls to the function)
    if(Trace){
      cat("\n We start the search of the value for coef=Info.max/If. \n \n")
    }        
    ## {{{ initialize key values to be updated in the following loop
    nwhile <- 0
    mycoefL0 <- mycoefL
    mycoefU <- mycoefMax
    mycoef <- mycoefU
    ## }}}
    ## {{{ Is the interval within which to search for coef large enough ?      
    if(Trace){
      cat("\n Check whether the interval within which to search for coef large enough. \n")
    }
    xx <- Method2_PC(rho_alpha=rho_alpha,
                       rho_beta=rho_beta,
                       alpha=alpha,
                       beta=beta,
                       Kmax=Kmax,         
                       Info.max=If*mycoefU,
                       InfoR.i=InfoR.i,
                       InfoR.d=InfoR.d,
                       delta=delta,
                       abseps=abseps,
                       toldiff=toldiff,
                       direction=direction,
                       Trace=FALSE,
                       sided=sided,
                       cMin=cMin)
    thediff <- abs(xx$boundaries[Kmax,"u.k"]-xx$boundaries[Kmax,"l.k"])
    ## }}}       
    if(thediff==0){
      ## {{{ if yes, we search coef
      mycoef <- (mycoefL + mycoefU)/2
      thediff <- 2*toldiff
      if(Trace){
        cat("\n we start the search within [",mycoefL,",",mycoefU,"] \n")
      }
      while(nwhile < nWhileMax & thediff>toldiff & abs(mycoefL-mycoefU)> tolcoef){
        nwhile <- nwhile + 1
        if(Trace){
          cat("\n Step :",nwhile,"(out of max.", nWhileMax,")")
        }
        xx <- Method2_PC(rho_alpha=rho_alpha,
                           rho_beta=rho_beta,
                           alpha=alpha,
                           beta=beta,
                           Kmax=Kmax,         
                           Info.max=If*mycoef,
                           InfoR.i=InfoR.i,
                           InfoR.d=InfoR.d,
                           delta=delta,
                           abseps=abseps,
                           toldiff=toldiff,
                           direction=direction,
                           sided=sided,
                           cMin=cMin)
        thediff <- abs(xx$boundaries[Kmax,"u.k"]-xx$boundaries[Kmax,"l.k"])
        if(thediff>toldiff){
          if(Trace){
            cat("\n Value coef=",mycoef,"is too small  \n")
            cat("\n coef=",mycoef,"leads to b.K-a.K=",thediff, "(whereas tol=",toldiff,") \n")
            ## cat("\n b.K=",xx$boundaries[K,"b.k"]," and a.K=",xx$boundaries[Kmax,"a.k"]," \n")
          }
          mycoefL <- (mycoefL+mycoefU)/2
          if(Trace){
            cat("\n we update the interval : [",mycoefL,",",mycoefU,"] \n")
          }
          mycoef <- (mycoefL+mycoefU)/2
        }
        if(thediff==0){
          if(Trace){
            cat("\n Value coef=",mycoef,"is too large  \n")
            cat("\n coef=",mycoef,"leads to b.K-a.K=",thediff, "(whereas tol=",toldiff,") \n")
            ## cat("\n b.K=",xx$boundaries[Kmax,"b.k"]," and a.K=",xx$boundaries[Kmax,"a.k"]," \n")
          }
          mycoefU <- (mycoefL+mycoefU)/2
          if(Trace){
            cat("\n we update the interval : [",mycoefL,",",mycoefU,"] \n")
          }
          mycoef <- (mycoefL+mycoefU)/2
          thediff <- 2*toldiff
        }
        if((thediff<=toldiff & thediff!=0) | abs(mycoefL-mycoefU)<= tolcoef ){
          if(Trace){
            cat("\n coef value FOUND : coef=",mycoef,"\n (leads to b.K-a.K=",thediff, " and tol.=",toldiff," and search interval length is=",abs(mycoefL-mycoefU),"and tol.=",tolcoef,")\n")
          }
          Info.max <- mycoef*If
          if(Trace){
            cat("\n Info.max computed as=",Info.max,"\n")
          }
          ## browser()
          ## print("Info.max created")
        }else{
          if(nwhile==nWhileMax){
            stop("Info.max could not be computed presicely enough : we need to allow for more iterations in the algorithm : you should probably call the function again with a larger value for nWhileMax.")
          }
        }
      }
      ## }}}
    }else{
      ## {{{ if no, we stop and explain why
      stop("The interval [mycoefL,mycoefMax]= [",mycoefL0,",",mycoefMax,"] is too small. You should probably call the function again with a larger value for mycoefMax and/or a lower (value >=1) for mycoefL \n")
      ## }}}
    }        
    
    ## }}}
  }else{
    mycoef <- Info.max/If
  }
  ## {{{ compute vector of means under the alternative H1
  
  #compute information at each analysis
  Info.i <- InfoR.i*Info.max
  Info.d <- InfoR.d*Info.max
  ImaxAnticipated <- Info.i>=Info.max | Info.d>=Info.max
  
  # compute the mean of the multivariate normal distribution under the alternative H1
  thetheta <- delta*sqrt(Info.i)
  ## }}} 
  ## {{{ case k=1 
  IncAlpha[1] <- ErrorSpend(I=Info.i[1],rho=rho_alpha,beta_or_alpha=alpha,Info.max=Info.max)
  IncBeta[1] <-  ErrorSpend(I=Info.i[1],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)
  
  uk[1] <- qnorm(p=1-IncAlpha[1],mean=0,sd=1)         # compute under the null (Ho)
  #thea[1] <- qnorm(p=IncBeta[1],mean=thetheta[1],sd=1)  # compute under the alternative (H1)

  #------------
  find.lk <- function(x){
    #calculate c corresponding to lk
    ck <- Method1(uk=uk[1],
                  lk=x,
                  Info.i=Info.i[1],
                  Info.d=Info.d[1],
                  Info.max=Info.max,
                  sided=sided,
                  cMin=cMin,
                  ImaxAnticipated=ImaxAnticipated[1],
                  rho_alpha=rho_alpha,
                  alpha=alpha,
                  bindingFutility = TRUE)
    
    #information matrix for first interim and decision analysis
    sigmaZk2 <- matrix(NA,ncol=2,nrow=2)
    sigmaZk2[1,1] <- sigmaZk[1,1]
    sigmaZk2[2,2] <- 1
    sigmaZk2[1,2] <- sigmaZk2[2,1] <- sqrt(Info.i[1]/Info.d[1])
    
    #probability to conclude futility
    pmvnorm(lower = c(uk[1],-Inf),
            upper = c(Inf,ck),
            mean=c(thetheta[1],delta*sqrt(Info.d[1])),
            sigma= sigmaZk2,
            abseps = abseps) +
    pmvnorm(lower = c(-Inf,-Inf),
            upper = c(x,ck),
            mean=c(thetheta[1],delta*sqrt(Info.d[1])),
            sigma= sigmaZk2,
            abseps = abseps) - IncBeta[1]
      
  }
  lk[1] <- uniroot(find.lk,lower=uk[1]-10,upper=uk[1])$root  #dirty solution to use -10 for lower bound
  ck[1] <- Method1(uk=uk[1],
                   lk=lk[1],
                   Info.i=Info.i[1],
                   Info.d=Info.d[1],
                   Info.max=Info.max,
                   sided=sided,
                   cMin=cMin,
                   ImaxAnticipated=ImaxAnticipated[1],
                   rho_alpha=rho_alpha,  #needs updating in Method1 to allow different rhos for beta and alpha
                   alpha=alpha,
                   bindingFutility = TRUE)
    #------------
  
  thealpha[1] <- IncAlpha[1]   
  thebeta[1] <- IncBeta[1]
  
  lk <- pmin(lk,uk) # just in case of over-running

  ## if(Trace){
  ## cat("\n a.1 computed as",lk[1],"and b.1 as",uk[1]," \n")
  ## }
  ## }}}
  ## {{{ loop over k >=2
  if(Kmax>1){
    for(k in 2:Kmax){
      if(!lk[k-1]==uk[k-1]){
        ## {{{ if  over-running has not occurred yet
        thealpha[k] <- ErrorSpend(I=Info.i[k],rho=rho_alpha,beta_or_alpha=alpha,Info.max=Info.max) 
        IncAlpha[k] <- thealpha[k] - thealpha[(k-1)]   
        thebeta[k] <- ErrorSpend(I=Info.i[k],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)  
        IncBeta[k] <- thebeta[k] - thebeta[(k-1)]   
        ## {{{ 
        ## {{{ u_k by solving what follows 
        uk[k] <- (uk[k-1] + lk[k-1])/2 # just to handle cases in which there is no root in what follows (when binding = TRUE )
        try(uk[k] <- uniroot(function(x){pmvnorm(lower = c(lk[1:(k-1)],x),
                                                 upper = c(uk[1:(k-1)],Inf),
                                                 mean=rep(0,k),
                                                 sigma= sigmaZk[1:k,1:k],
                                                 abseps = abseps) - IncAlpha[k]},
                                 lower = lk[k-1],
                                 upper = uk[k-1],
                                 tol = abseps)$root, silent = TRUE)
        IsbkOK <- !(uk[k]==((uk[k-1] + lk[k-1])/2))
        ## }}}
        ## {{{ a_k by solving what follows 
        if(IsbkOK){
          if(k!=Kmax){
            find.lkk <- function(x){
              #calculate c corresponding to lk
              ck <- Method1(uk=uk[1:k],
                            lk=c(lk[1:(k-1)],x),
                            Info.i=Info.i[1:k],
                            Info.d=Info.d[k],
                            Info.max=Info.max,
                            sided=sided,
                            cMin=cMin,
                            ImaxAnticipated=ImaxAnticipated[k],
                            rho_alpha=rho_alpha,
                            alpha=alpha,
                            bindingFutility = TRUE)
                
              #information matrix for first interim and decision analysis
              sigmaZk2 <- matrix(NA,ncol=k+1,nrow=k+1)
              sigmaZk2[1:k,1:k] <- sigmaZk[1:k,1:k]
              sigmaZk2[k+1,k+1] <- 1
              sigmaZk2[1:k,k+1] <- sigmaZk2[k+1,1:k] <- sqrt(Info.i[1:k]/Info.d[k])
                
              #probability to conclude futility
              pmvnorm(lower = c(lk[1:(k-1)],uk[k],-Inf),
                      upper = c(uk[1:(k-1)],Inf,ck),
                      mean=c(thetheta[1:k],delta*sqrt(Info.d[k])),
                      sigma= sigmaZk2,
                      abseps = abseps) +
              pmvnorm(lower = c(lk[1:(k-1)],-Inf,-Inf),
                      upper = c(uk[1:(k-1)],x,ck),
                      mean=c(thetheta[1:k],delta*sqrt(Info.d[k])),
                      sigma= sigmaZk2,
                      abseps = abseps) - IncBeta[k]
              
            }
            lk[k] <- try(uniroot(find.lkk,lower=uk[k]-10,upper=uk[k])$root)

            if(!inherits(lk[k], "try-error")){
              ck[k] <- Method1(uk=uk[1:k],
                             lk=lk[1:k],
                             Info.i=Info.i[1:k],
                             Info.d=Info.d[k],
                             Info.max=Info.max,
                             sided=sided,
                             cMin=cMin,
                             ImaxAnticipated=ImaxAnticipated[k],
                             rho_alpha=rho_alpha,
                             alpha=alpha,
                             bindingFutility = TRUE)
            } else {
              lk[k] <- uk[k] # just to handle cases in which there is no root
              if(inherits(lk[k],"try-error")){warning(paste0("try-error for calculation of lk[",k,"]"))}
            }
          }
          if(k==Kmax){
            lk[k] <- uk[k] # just to handle cases in which there is no root
            
            try(lk[k] <- uniroot(function(x){pmvnorm(lower = c(lk[1:(k-1)],-Inf),
                                                     upper = c(uk[1:(k-1)],x),
                                                     mean=thetheta[1:k],
                                                     sigma= sigmaZk[1:k,1:k],
                                                     abseps = abseps) - IncBeta[k]},
                                   lower = lk[k-1], 
                                   upper = uk[k], 
                                   tol = abseps)$root, silent = TRUE)
            if(inherits(lk[k],"try-error")){warning("try-error for calculation of lk[Kmax]")}
                
          }
        }else{
            lk[k] <- (uk[k-1] + lk[k-1])/2 # just to handle cases in which there is no root in what is above
        }
          ## }}}
        ## {{{ to deal with over-running (see chapter Jennison) 
        lk <- pmin(lk,uk)
        ## }}}
      }else{
        ## {{{ if  over-running has already occurred
        lk[k:Kmax] <- lk[k-1]
        uk[k:Kmax] <- uk[k-1]
        ## }}}
      }
    }
  }
  ## }}}
  
  if(direction=="greater"){
    lk <- -lk
    uk <- -uk
    delta <- -delta
  }
  
  ## {{{ create  output
  d <- data.frame(l.k=lk,
                  u.k=uk,
                  c.k=ck,
                  Type.I.Error=thealpha,
                  Type.II.Error=thebeta,
                  Inc.Type.I=IncAlpha,
                  Inc.Type.II=IncBeta,
                  I.k=Info.i
  )
  out <- list(boundaries=d,
              rho_alpha=rho_alpha,
              rho_beta=rho_beta,
              alpha=alpha,
              beta=beta,
              Kmax=Kmax,
              If=If,
              Info.max=Info.max,
              Info.i=Info.i,
              Info.d=Info.d,
              delta=delta,
              coef=mycoef,
              abseps=abseps,
              toldiff=toldiff,
              direction=direction,
              sided=sided,
              binding=TRUE,
              cMin=cMin)
  class(out) <- "delayedGSD"
  ## }}}
  .Random.seed <<- old # restore the current seed (before the call to the function)
  out
}
