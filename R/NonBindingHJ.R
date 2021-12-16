#updated function from previous tutorial such that it uses same naming and errorspending functions
#updated function such that it can handle any specified information levels rather than just equally spaced

#next step: implement HJ non-binding rule
#TBD: do we want to use our own function instead of gsDesign?

#' @title Calculate boundaries for a group sequential design
#' @description Calculate boundaries for a group sequential design without delayed endpoints based on planned and/or observed information using an error spending approach. The function can also give the boundaries for a non-binding futility rule for delayed endpoints following a method proposed by Hampson and Jennison.
#' 
#' 
#' 
#' @param rho_alpha rho parameter of the rho-family spending functions (Kim-DeMets) for alpha
#' @param rho_beta rho parameter of the rho-family spending functions (Kim-DeMets) for beta
#' @param alpha type I error
#' @param beta type II error
#' @param kMax max number of analyses (including final)
#' @param Info.max maximum information needed for given beta (type II-error), theta (expected difference), alpha (type I-error) and Kmax. It can be given if it is known. Otherwise it is computed from the  values given for alpha, beta, theta and Kmax.
#' @param Info.i Expected or observed (wherever possible) information at the interim and final analyses 1:Kmax
#' @param Info.d Expected or observed information at all potential decision analyses 1:(Kmax-1)
#' @param theta expected effect under the alternative (should be on the scale of the test statistc for which If and Info.max relate to one over the variance, e.g. theta=expected log(Hazard ratio))
#' @param abseps tolerance for precision when finding roots or computing integrals
#' @param binding binding or non-binding futility boundary (i.e FALSE if it is allowed to continue after crossing the futility boundary)
#' @param binding_type type of non-binding futility rule to use. "HJ" corresponds to the method proposed by Hampson and Jennison, "BBO" corresponds to the method proposed by Baayen, Blanche and Ozenne
#' @param direction greater is for Ho= theta > 0, "smaller" is for Ho= theta < 0 (note that in Jennison and Turnbull's book chapter (2013) they consider smaller)
#' @param Trace Used only if Info.max=NULL. Whether to print informations to follow the progression of the (root finding) algorithm to compute Info.max (from  alpha, beta, theta and Kmax).
#' @param nWhileMax Used only if Info.max=NULL. Maximum number of steps in the (root finding) algorithm to compute Info.max (from  alpha, beta, theta and Kmax)
#' @param toldiff Used only if Info.max=NULL. Maximum tolerated difference between lower and upper bounds at anaylis Kmax (which souhld be zero), in the root finding algorithm, to find the value of Info.max
#' @param tolcoef Used only if Info.max=NULL. Maximum tolerated difference before stopping the search (in the root finding algorithm), between two successive values for the multiplier coeficient 'coef' such that Info.max=coef*If  (some values for coef are given in Table 7.6 page 164 Jennison's book. The value of "If" (which stands for Information for Fixed design) corresponds to the information we need if Kmax=1)
#' @param mycoefMax Used only if Info.max=NULL. Upper limit of the interval of values in which we search for the multiplier coeficient 'coef' such that Info.max=coef*If (in the root finding algorithm).
#' @param mycoefL Used only if Info.max=NULL. Lower limit of the interval (see mycoefMax)
#' @param myseed seed for producing reproducible results. Because we call functions which are based on Monte-Carlo compuation (pmvnorm)
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
#'                     gammaA=2,  #rho parameter for alpha error spending function
#'                     gammaB=2,  #rho parameter for beta error spending function
#'                     method=1,  #use method 1 or 2 from paper H&J
#'                     delta=1.5,  #effect that the study is powered for
#'                     InfoR.d=0.65,
#'                     bindingFutility=FALSE)
#'
#' b12 <- NonBindingHJ(Kmax=2,Info.max=b1$Info.max,theta=1.5,binding=FALSE,alpha=0.025,Info.i=c(0.6,1)*b1$Info.max)
#' 
#' all.equal(b1$uk, b12$boundaries[,"b.k"])
#' all.equal(b1$lk, b12$boundaries[,"a.k"])

## * NonBindingHJ (code)
NonBindingHJ <- function(rho_alpha=2,          # rho parameter of the rho-family spending functions (Kim-DeMets) for alpha
                         rho_beta=2,           # rho parameter of the rho-family spending functions (Kim-DeMets) for beta
                         alpha=0.025,          # Type-I error (overall)
                         beta=0.2,             # Type-II error (overall)
                         Kmax,                 # number of planned analyses (including the final analysis)
                         Info.max=NULL,        # Info.max, i.e. maximum information needed for given beta (type II-error), theta (expected difference), alpha (type I-error) and Kmax. It can be given if it is known. Otherwise it is computed from the  values given for alpha, beta, theta and Kmax.
                         Info.i=NULL,          # Expected or observed (wherever possible) information at the interim and final analyses 1:Kmax
                         Info.d=NULL,          # Expected or observed information at all potential decision analyses 1:(Kmax-1)
                         theta=0,              # expected effect under the alternative (should be on the scale of the test statistc for which If and Info.max relate to one over the variance, e.g. theta=expected log(Hazard ratio))
                         abseps = 1e-06,       # tolerance for precision when finding roots or computing integrals
                         binding=TRUE,         # binding or non-binding futility boundary (i.e FALSE if it is allowed to continue after crossing the futility boundary)
                         binding_type="BBO",   # type of non-binding futility rule to use. "HJ" corresponds to the method proposed by Hampson and Jennison, "BBO" corresponds to the method proposed by Baayen, Blanche and Ozenne
                         direction="smaller",  # greater is for Ho= theta > 0, "smaller" is for Ho= theta < 0 (note that in Jennison and Turnbull's book chapter (2013) they consider smaller)
                         Trace=FALSE,          # Used only if Info.max=NULL. Whether to print informations to follow the progression of the (root finding) algorithm to compute Info.max (from  alpha, beta, theta and Kmax).
                         nWhileMax=30,         # Used only if Info.max=NULL. Maximum number of steps in the (root finding) algorithm to compute Info.max (from  alpha, beta, theta and Kmax)
                         toldiff= 1e-05,       # Used only if Info.max=NULL. Maximum tolerated difference between lower and upper bounds at anaylis Kmax (which souhld be zero), in the root finding algorithm, to find the value of Info.max
                         tolcoef= 1e-04,       # Used only if Info.max=NULL. Maximum tolerated difference before stopping the search (in the root finding algorithm), between two successive values for the multiplier coeficient 'coef' such that Info.max=coef*If  (some values for coef are given in Table 7.6 page 164 Jennison's book. The value of "If" (which stands for Information for Fixed design) corresponds to the information we need if Kmax=1)
                         mycoefMax= 1.2,       # Used only if Info.max=NULL. Upper limit of the interval of values in which we search for the multiplier coeficient 'coef' such that Info.max=coef*If (in the root finding algorithm).
                         mycoefL=1,            # Used only if Info.max=NULL. Lower limit of the interval (see mycoefMax)
                         myseed=2902           # seed for producing reproducible results. Because we call functions which are based on Monte-Carlo compuation (pmvnorm)
){

  ## {{{ set seed
  old <- .Random.seed # to save the current seed
  set.seed(myseed)
  ## }}}
  ## {{{ preliminaries
  mycoef <- NULL # initialize as needed for output
  if(direction=="smaller"){
    thea <- rep(-Inf,Kmax)
    theb <- rep(Inf,Kmax)     
  }else{
    if(direction=="greater"){
      thea <- rep(Inf,Kmax)
      theb <- rep(-Inf,Kmax)
    }else{
      stop("direction should be either greater or smaller")}
  }

  if( (direction=="smaller" & theta<0) | (direction=="greater" & theta>0)){
    stop("The values given for arguments direction and theta are inconsistent.\n When direction=smaller, theta should be positive.\n When direction=greater, theta should be negative.")
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
      sigmaZk[i,j] <- sqrt(Info.i[i]/Info.i[j])
      sigmaZk[j,i] <- sqrt(Info.i[i]/Info.i[j])
    }
  }
  # compute If (see Jennison book page 87
  If <- (qnorm(1-alpha)+qnorm(1-beta))^2/theta^2
  if(Trace){
    cat("\n If computed as =",If,"\n")
  }
  ## }}}
  if(is.null(Info.max)){
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
    xx <- NonBindingHJ(rho_alpha=rho_alpha,
                         rho_beta=rho_beta,
                         alpha=alpha,
                         beta=beta,
                         Kmax=Kmax,         
                         Info.max=If*mycoefU,
                         Info.i=Info.i,
                         Info.d=Info.d,
                         theta=theta,
                         abseps=abseps,
                         toldiff=toldiff,
                         binding=binding,
                         binding_type=binding_type,
                         direction=direction,
                         Trace=FALSE)
    thediff <- abs(xx$boundaries[Kmax,"b.k"]-xx$boundaries[Kmax,"a.k"])
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
        xx <- NonBindingHJ(rho_alpha=rho_alpha,
                             rho_beta=rho_beta,
                             alpha=alpha,
                             beta=beta,
                             Kmax=Kmax,         
                             Info.max=If*mycoef,
                             Info.i=Info.i,
                             Info.d=Info.d,
                             theta=theta,
                             abseps=abseps,
                             toldiff=toldiff,
                             binding=binding,
                             binding_type=binding_type,
                             direction=direction)
        thediff <- abs(xx$boundaries[Kmax,"b.k"]-xx$boundaries[Kmax,"a.k"])
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
  # compute the mean of the multivariate normal distribution under the alternative H1
  thetheta <- theta*sqrt(Info.i)
  ## }}} 
  ## {{{ case k=1 
  IncAlpha[1] <- ErrorSpend(I=Info.i[1],rho=rho_alpha,beta_or_alpha=alpha,Info.max=Info.max)
  IncBeta[1] <-  ErrorSpend(I=Info.i[1],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)
  if(direction=="smaller"){
    theb[1] <- qnorm(p=1-IncAlpha[1],mean=0,sd=1)         # compute under the null (Ho)
    thea[1] <- qnorm(p=IncBeta[1],mean=thetheta[1],sd=1)  # compute under the alternative (H1)
  }else{
    if(direction=="greater"){
      theb[1] <- qnorm(p=IncAlpha[1],mean=0,sd=1)         # compute under the null (Ho)
      thea[1] <- qnorm(p=1-IncBeta[1],mean=thetheta[1],sd=1)  # compute under the alternative (H1)
    }else{
      stop("direction should be either greater or smaller")
    }
  }
  
  thealpha[1] <- IncAlpha[1]   
  thebeta[1] <- IncBeta[1]
  if(direction=="smaller"){
    thea <- pmin(thea,theb) # just in case of over-running
  }else{
    theb <- pmin(thea,theb) # just in case of over-running
  }
  ## if(Trace){
  ## cat("\n a.1 computed as",thea[1],"and b.1 as",theb[1]," \n")
  ## }
  ## }}}
  ## {{{ loop over k >=2
  if(Kmax>1){
    for(k in 2:Kmax){
      if(!thea[k-1]==theb[k-1]){
        ## {{{ if  over-running has not occurred yet
        thealpha[k] <- ErrorSpend(I=Info.i[k],rho=rho_alpha,beta_or_alpha=alpha,Info.max=Info.max) 
        IncAlpha[k] <- thealpha[k] - thealpha[(k-1)]   
        thebeta[k] <- ErrorSpend(I=Info.i[k],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)  
        IncBeta[k] <- thebeta[k] - thebeta[(k-1)]   
        if(binding){
          ## {{{ 
          if(direction=="smaller"){
            ## {{{ b_k by solving what follows 
            theb[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what follows (when binding = TRUE )
            try(theb[k] <- uniroot(function(x){pmvnorm(lower = c(thea[1:(k-1)],x),
                                                       upper = c(theb[1:(k-1)],Inf),
                                                       mean=rep(0,k),
                                                       sigma= sigmaZk[1:k,1:k],
                                                       abseps = abseps) - IncAlpha[k]},
                                   lower = thea[k-1],
                                   upper = theb[k-1],
                                   tol = abseps)$root, silent = TRUE)
            IsbkOK <- !(theb[k]==((theb[k-1] + thea[k-1])/2))
            ## }}}
            ## {{{ a_k by solving what follows 
            if(IsbkOK){
              thea[k] <- theb[k] # just to handle cases in which there is no root in what follows 
              try(thea[k] <- uniroot(function(x){pmvnorm(lower = c(thea[1:(k-1)],-Inf),
                                                         upper = c(theb[1:(k-1)],x),
                                                         mean=thetheta[1:k],
                                                         sigma= sigmaZk[1:k,1:k],
                                                         abseps = abseps) - IncBeta[k]},
                                     lower = thea[k-1], 
                                     upper = theb[k], 
                                     tol = abseps)$root, silent = TRUE)
            }else{
              thea[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what is above
            }
            ## }}}
          }
          if(direction=="greater"){
            ## {{{ b_k by solving what follows 
            theb[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what follows  (when binding = TRUE )
            try(theb[k] <- uniroot(function(x){pmvnorm(lower = c(theb[1:(k-1)],-Inf),
                                                       upper = c(thea[1:(k-1)],x),
                                                       mean=rep(0,k), 
                                                       sigma= sigmaZk[1:k,1:k],
                                                       abseps = abseps) - IncAlpha[k]},
                                   lower = theb[k-1],
                                   upper = thea[k-1],
                                   tol = abseps)$root, silent = TRUE)
            IsbkOK <- !(theb[k]==((theb[k-1] + thea[k-1])/2))
            ## }}}
            ## {{{ a_k by solving what follows 
            if(IsbkOK){
              thea[k] <- theb[k] # just to handle cases in which there is no root in what follows                            
              try(thea[k] <- uniroot(function(x){pmvnorm(lower = c(theb[1:(k-1)],x),
                                                         upper = c(thea[1:(k-1)],Inf),
                                                         mean=thetheta[1:k],
                                                         sigma= sigmaZk[1:k,1:k],
                                                         abseps = abseps) - IncBeta[k]},
                                     lower = theb[k],
                                     upper = thea[k-1],
                                     tol = abseps)$root, silent = TRUE)
            }else{
              thea[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what is above
            }
            ## }}}
          }
          ## }}}
        }else if(binding_type=="BBO"){
            ## {{{ 
            if(direction=="smaller"){
              ## {{{ b_k by solving what follows
              theb[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what follows
              try(theb[k] <- uniroot(function(x){pmvnorm(lower = c(rep(-Inf, k-1),x),
                                                         upper = c(theb[1:(k-1)],Inf),                                                                   
                                                         mean=rep(0,k),
                                                         sigma= sigmaZk[1:k,1:k],
                                                         abseps = abseps) - IncAlpha[k]},
                                     lower = thea[1],
                                     upper = theb[1],
                                     tol = abseps)$root, silent = TRUE)
              ## browser()
              IsbkOK <- !(theb[k]==((theb[k-1] + thea[k-1])/2))
              ## }}}
              ## {{{ a_k by solving what follows
              if(IsbkOK){
                thea[k] <- theb[k] # just to handle cases in which there is no root in what follows 
                try(thea[k] <- uniroot(function(x){pmvnorm(lower = c(thea[1:(k-1)],-Inf),
                                                           upper = c(theb[1:(k-1)],x),
                                                           mean=thetheta[1:k],
                                                           sigma= sigmaZk[1:k,1:k],
                                                           abseps = abseps) - IncBeta[k]},
                                       lower = thea[1],
                                       upper = theb[1],
                                       tol = abseps)$root, silent = TRUE)
              }else{
                thea[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what is above
              }
              ## }}}
            }
            if(direction=="greater"){
              ## {{{ b_k by solving what follows
              ## browser()
              theb[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what follows
              try(theb[k] <- uniroot(function(x){pmvnorm(lower = c(theb[1:(k-1)],-Inf),
                                                         upper = c(rep(Inf, k-1),x),
                                                         mean=rep(0,k), #thetheta[1:k],
                                                         sigma= sigmaZk[1:k,1:k],
                                                         abseps = abseps) - IncAlpha[k]},
                                     lower = theb[k-1],
                                     upper = thea[k-1],
                                     tol = abseps)$root, silent = TRUE)
              IsbkOK <- !(theb[k]==((theb[k-1] + thea[k-1])/2))
              ## }}}
              ## {{{ a_k by solving what follows
              if(IsbkOK){
                thea[k] <- theb[k] # just to handle cases in which there is no root in what follows                            
                try(thea[k] <- uniroot(function(x){pmvnorm(lower = c(theb[1:(k-1)],x),
                                                           upper = c(thea[1:(k-1)],Inf),
                                                           mean=thetheta[1:k],
                                                           sigma= sigmaZk[1:k,1:k],
                                                           abseps = abseps) - IncBeta[k]},
                                       lower = theb[1],
                                       upper = thea[1],
                                       tol = abseps)$root, silent = TRUE)
              }else{
                thea[k] <- (theb[k-1] + thea[k-1])/2 # just to handle cases in which there is no root in what is above
              }
              ## }}}
            }
            ## }}}
          } else if(binding_type=="HJ"){
            stop("not implemented yet")
          } else {
            stop("Please specify binding_type is HJ or BBO")
          }
        ## {{{ to deal with over-running (see chapter Jennison) 
        if(direction=="smaller"){
          thea <- pmin(thea,theb)
        }else{
          theb <- pmin(thea,theb)
        }
        ## }}}
        ## }}}
      }else{
        ## {{{ if  over-running has already occurred
        thea[k:Kmax] <- thea[k-1]
        theb[k:Kmax] <- theb[k-1]
        ## }}}
      }
    }
  }
  ## }}}
  ## {{{ create  output
  d <- data.frame(a.k=thea,
                  b.k=theb,
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
              theta=theta,
              coef=mycoef,
              abseps=abseps,
              toldiff=toldiff,
              binding=binding,
              binding_type=binding_type,
              direction=direction)
  class(out) <- "SeqCR"
  ## }}}
  .Random.seed <<- old # restore the current seed (before the call to the function)
  out
}
