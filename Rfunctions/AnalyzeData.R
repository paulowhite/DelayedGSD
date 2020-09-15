### AnalyzeData.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 25 2020 (09:23) 
## Version: 
## Last-Updated: Aug 28 2020 (17:14) 
##           By: Paul Blanche
##     Update #: 32
#----------------------------------------------------------------------
## 
### Commentary: 
##
## Analyze available data: provide estimate, s.e., Information nd p-value.
##
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

AnalyzeData <- function(d){

    require(nlme)
    N.fw <- length(grep("t",names(d)))
    # transform to long format
    d$baseline <- d$X1
    d$X1 <- NULL
    d$t1 <- NULL
    
    long <- reshape(d, direction='long', varying=paste0("X",2:N.fw), idvar='id', v.names='X')
    long$id <- factor(long$id)
    long$Z <- factor(long$Z)
    long$time.factor <- factor(long$time)
    long$time.factor <- relevel(long$time.factor,ref=N.fw-1)
    head(long[order(long$id),],n=10)
    summary(long)
    
    ctrl <- glsControl(opt='optim')
    
    m <- gls(X~baseline*time.factor + Z*time.factor,
             data=long,
             correlation=corSymm(form=~1|id), 
             varIdent(form=~1|time),method="REML",
             na.action = na.exclude)

    ## ctrl <- glsControl(opt='optim')
    ## fit0 <- gls(X~baseline*time.factor + Z*time.factor,
    ## data=long, 
    ## correlation=corSymm(form=~time|id), 
    ## weights=varIdent(form=~1|time), 
    ## na.action = na.exclude, control=ctrl)
    #design matrix

    ## Xmat <- rbind(c(1,mean(d$baseline),0,1,0,0,baseline_mean,0,0,0,0,0), #active trt, w16
                  ## c(1,baseline_mean,1,0,0,1,baseline_mean,0,0,1,0,0)) #control trt, w16
    
    ## cmat <- Xmat%*%m$varBeta%*%t(Xmat)   #covariance matrix predicted lsmeans   
    ## means <- coef(m)%*%t(Xmat)   #lsmeans
    ## meandiff <- means[1]-means [2]   #lsmeans difference
    ## SEmeandiff <- sqrt(cmat [1,1]+cmat [2,2]-2*cmat [1,2])   #SE of mean difference
    ## pval_2_sided <- 2*(pnorm(-abs(meandiff/SEmeandiff)))  #2-sided p-value
    ## pval_1_sided <- pnorm(meandiff/SEmeandiff, lower.tail=TRUE)    #1-sided p-value

   
    out <- list(d.long=long,
                fit=m,
                pavlue=2*(pnorm(-abs(coef(m)["Z1"]/sqrt(m$varBeta["Z1","Z1"])))),
                estimate=coef(m)["Z1"],
                se=sqrt(m$varBeta["Z1","Z1"]),
                Info=1/m$varBeta["Z1","Z1"])
    out
}

#----------------------------------------------------------------------
### AnalyzeData.R ends here
