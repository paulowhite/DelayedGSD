### AnalyzeData.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 25 2020 (09:23) 
## Version: 
## Last-Updated: Mar  9 2021 (12:54) 
##           By: Paul Blanche
##     Update #: 44
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

#' @param d dataset
#' @param ddf [character] method used to compute the degrees of freedom of the Wald statistic.
#' Can be \code{satterthwaite} or \code{nlme}.
AnalyzeData <- function(d, ddf = "nlme"){

    require(nlme)

    ddf <- match.arg(ddf, choices = c("nlme","satterthwaite"))
    N.fw <- length(grep("t",names(d)))
    # transform to long format
    d$baseline <- d$X1
    d$X1 <- NULL
    d$t1 <- NULL
    
    long <- reshape(d, direction='long', varying=paste0("X",2:N.fw), idvar='id', v.names='X')
    long$id <- factor(long$id)
    long$Z <- factor(long$Z)
    long$time.factor <- factor(long$time)
    long$time.factor <- relevel(long$time.factor,ref=N.fw-1) # relvel to make the treatment coefficient have the interpretation of threatment effect at end of follow-up
    ## head(long[order(long$id),],n=10)
    ## summary(long)
    
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


    if(ddf=="nlme"){
        out <- list(d.long=long,
                    fit=m,
                    pvalue=2*(pnorm(-abs(coef(m)["Z1"]/sqrt(m$varBeta["Z1","Z1"])))),
                    estimate=coef(m)["Z1"],
                    se=sqrt(m$varBeta["Z1","Z1"]),
                    Info=1/m$varBeta["Z1","Z1"])
    }else{
        require(emmeans)
        groupTest <- emmeans::emmeans(m, specs = ~Z|time.factor, data = long[!is.na(long$X),,drop=FALSE])
        e.satterthwaite <- summary(pairs(groupTest, reverse = TRUE), by = NULL, infer = TRUE, adjust = "none")
        index <- which(e.satterthwaite$time.factor==levels(long$time.factor)[1])
        
        out <- list(d.long=long,
                    fit=m,
                    pvalue=e.satterthwaite$p.value[index],
                    estimate=e.satterthwaite$estimate[index],
                    se=e.satterthwaite$SE[index],
                    Info=1/e.satterthwaite$SE[index]^2)

        if(abs(out$estimate-coef(m)["Z1"])>1e-12){
            warning("Something went wrong when extracting the coefficient from emmeans \n")
        }
    }
    
    return(out)
}

#----------------------------------------------------------------------
### AnalyzeData.R ends here
