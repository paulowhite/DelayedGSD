### GenData.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 25 2020 (08:14) 
## Version: 
## Last-Updated: Aug 28 2020 (16:32) 
##           By: Paul Blanche
##     Update #: 59
#----------------------------------------------------------------------
## 
### Commentary: 
##  
##  Generate longitudnal (full) data.
##  
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


##' @title Generate longitudnal (full) data
##' @description yy
##' @param n sample size
##' @param N.fw  number of measurements, i.e. 1 (for baseline) + number of follow-up measurements (equally spaced)
##' @param rand.block for block randomization
##' @param sd1 sd of primary outcome at end of follow-up (main outcome)
##' @param delta1 treatment effect on primary outcome (main outcome)
##' @param ar accrual rate (unit is per time between the equally spaced visits)
##' @param cor.12.1 correlation between outcome at baseline and at first visit  (main outcome)
##' @param cor.ij.1 correlation between outcome at two consecutive measurements  (main outcome)
##' @param initial.sd1 sd of primary outcome at baseline (main outcome). We suppose it evoluates linearly to sd1.
##' @param seed
##' @details zz
##' @return ff
##' @author Paul Blanche
##' 
##' @examples
##' x <- GenData()
##' head(x$d,n=20)
##' @export
GenData <- function(n=50, 
                    N.fw=4,
                    rand.block=c(1,1,0,0),
                    sd1=3,
                    delta1=1.5,
                    ar=10, 
                    cor.12.1=-0.3,
                    cor.ij.1=0.34,
                    initial.sd1=2.2,
                    seed=24082020
                    ){
    #--
    require(mvtnorm)
    #--
    set.seed(seed)   
    d <- data.frame(id=1:n, # id
                    Z=blockrand(seed,rand.block,n))
                    ## Z=c(rep(1,floor(n/2)),rep(0,n-floor(n/2)))[sample(1:n)])     # treatment
    for(k in 1:N.fw){
        d[[paste0("X",k)]] <- rep(NA,n) # initialize measurement values at each follow-up
        ## d[[paste0("t",k)]] <- rep(NA,n) # initialize time of each follow-up measurement
    }
    # correlation matrix of outcome at baseline and each follow-up time (can be made an independent subfunction)
    cor.matrix <- matrix(NA,ncol=N.fw,nrow=N.fw)
    diag(cor.matrix) <- 1
    cor.matrix[1,2] <- cor.matrix[2,1] <- cor.12.1
    if(N.fw>2){
        for(k in 3:N.fw){
            for(l in 1:(k-1)){
                cor.matrix[k,l] <- cor.matrix[l,k] <- cor.ij.1^(k-2)
            }
        }
        cor.matrix[1,3:N.fw] <- cor.matrix[3:N.fw,1] <- cor.matrix[1,3:N.fw]*cor.12.1
    }
    cor.matrix
    # make sd vector
    allsd <- initial.sd1 + ((1:N.fw)-1)*(sd1-initial.sd1)/(N.fw-1)
    # Combine cor matrix and var vector to var-cov matrix
    cor.matrix%*%diag(allsd)
    diag(allsd)%*%cor.matrix%*%diag(allsd)
    vcovmat <- diag(allsd)%*%cor.matrix%*%diag(allsd)
    # Generate outcome at baseline and each follow-up visit
    allmeans0 <- rep(0,N.fw)
    ## allmeans1 <- 0 + ((1:N.fw)-1)*(delta1-0)/(N.fw-1) # assume linear trend
    allmeans1 <- rep(delta1,N.fw)  # assume same effect at all follow-up times
    d[d$Z==0,paste0("X",1:N.fw)] <- rmvnorm(n=sum(d$Z==0),mean=allmeans0,sigma=vcovmat)
    d[d$Z==1,paste0("X",1:N.fw)] <- rmvnorm(n=sum(d$Z==1),mean=allmeans1,sigma=vcovmat)    
    # Generate time at inclusion (time since study start) and deduce time at each follow-up
    ## d[,"t1"] <- (1:n)/ar
    for(k in 1:N.fw){
        d[[paste0("t",k)]] <- (1:n)/ar + k-1 # initialize time of each follow-up measurement
    }    
    ## head(d)
    list(cor=cor.matrix,sd=allsd,vcov=vcovmat,mean1=allmeans1,mean0=allmeans0,d=d)
}

#----------------------------------------------------------------------
### GenData.R ends here
