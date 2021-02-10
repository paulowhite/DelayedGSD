### FCT_fixedSimpleGSD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  3 2021 (15:10) 
## Version: 
## Last-Updated: feb 10 2021 (12:35) 
##           By: Brice Ozenne
##     Update #: 27
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * simData
#' @examples
#' dt.tempo <- simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.9, n.batch = 4)
simData <- function(n, sigma2, mu0, mu1, rho, n.batch){
    require(data.table)
    require(mvtnorm)
    
    Sigma <-  (diag(1-rho,2,2) + rho)*sigma2
    
    dt <- NULL

    for(iBatch in 1:n.batch){ ## iBatch <- 1

        iValue.control <- mvtnorm::rmvnorm(n = n, mean = rep(mu0,2), sigma = Sigma)
        iValue.treatment <- mvtnorm::rmvnorm(n = n, mean = rep(mu1,2), sigma = Sigma)
        
        for(iTime in 0:n.batch){
            if(iTime == (iBatch-1)){
                dt <- rbind(dt,
                            data.table(group = "control",
                                       id = paste0("c",(iBatch-1)*n+1:n),
                                       batch = iBatch,
                                       time = iTime,
                                       proxy = iValue.control[,1],
                                       outcome = as.numeric(NA)),
                            data.table(group = "treatment",
                                       id = paste0("t",(iBatch-1)*n+1:n),
                                       batch = iBatch,
                                       time = iTime,
                                       proxy = iValue.treatment[,1],
                                       outcome = as.numeric(NA))
                            )
            }else if(iTime > (iBatch-1)){
                dt <- rbind(dt,
                            data.table(group = "control",
                                       id = paste0("c",(iBatch-1)*n+1:n),
                                       batch = iBatch,
                                       time = iTime,
                                       proxy = iValue.control[,1],
                                       outcome = iValue.control[,2]),
                            data.table(group = "treatment",
                                       id = paste0("t",(iBatch-1)*n+1:n),
                                       batch = iBatch,
                                       time = iTime,
                                       proxy = iValue.treatment[,1],
                                       outcome = iValue.treatment[,2])
                            )
            }else{
                dt <- rbind(dt,
                            data.table(group = "control",
                                       id = paste0("c",(iBatch-1)*n+1:n),
                                       batch = iBatch,
                                       time = iTime,
                                       proxy = NA,
                                       outcome = NA),
                            data.table(group = "treatment",
                                       id = paste0("t",(iBatch-1)*n+1:n),
                                       batch = iBatch,
                                       time = iTime,
                                       proxy = NA,
                                       outcome = NA)
                            )
            }
        }
        
    }
    dt[,group := relevel(as.factor(group),"control")]
    attr(dt,"param") <- c(n=n, sigma2=sigma2, mu0=mu0, mu1=mu1, rho=rho, n.batch=n.batch)

    return(dt)
}

## * analyzeData
#' @examples
#' set.seed(10)
#' analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.5, n.batch = 4))$info
#' analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0, n.batch = 4))$info
#' analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.999, n.batch = 4))$info
analyzeData <- function(data, fast = TRUE){
    require(nlme)
    require(data.table)
    source("./Rfunctions/ttest.R")
    source("./Rfunctions/getInformation.R")

    out <- list()
    ## *** prepare
    n.obs.outcome <- data[!is.na(outcome),.N,by="time"]
    times <- unique(data$time)
    interim.time <- 2
    decision.time <- 3
    final.time <- 4

    if(any(c(interim.time, final.time) %in% times == FALSE)){stop("Incorrect time specification \n")}

    ## *** interim analysis (2 months, stop or continue)
    data.interim <- data[time == interim.time & (batch-1) <= interim.time]
    n.decision <- n.obs.outcome[time==decision.time,N]

    ## **** t-test
    out$interim.ttest <- ttest(outcome ~ group, data = data.interim)
    if(fast){
        out$info <- c(out$info,"info.ttest" = 1/out$interim.ttest$stderr^2)
        out$info <- c(out$info,"info.ttest2" = 1/out$interim.ttest$stderr^2 * n.decision/n.obs.outcome[time==interim.time,N])
    }else{
        out$info <- c(out$info,"info.ttest" = getInformation(out$interim.ttest, type = "estimation"))
        out$info <- c(out$info,"info.ttest2" = getInformation(out$interim.ttest, type = "prediction")) ## as.double(out$info["info.ttest"]) * pc.interim
    }
    
    ## **** mixed model
    dataL.interim <- melt(data.interim, id.vars = c("group","id","batch","time"), value.name = "value", variable.name = "type")
    dataL.interim$type <- relevel(dataL.interim$type, "outcome")
    out$interim.lmm <- gls(value ~ type*group,
                           data = dataL.interim,
                           weights = varIdent(form=~1|type*group),
                           correlation = corCompSymm(form=~1|id),
                           na.action = na.omit
                           )

    
    
    rho <- as.double(coef(out$interim.lmm$modelStruct$corStruct, unconstrain = FALSE))
    sigma <- (coef(out$interim.lmm$modelStruct$varStruct, unconstrain = FALSE, allCoef = TRUE)*sigma(out$interim.lmm))^2
    
    if(fast){
        out$info <- c(out$info,"info.lmm" = 1/vcov(out$interim.lmm)["grouptreatment","grouptreatment"])
    }else{
        out$info <- c(out$info,"info.lmm" = getInformation(out$interim.lmm, data = dataL.interim, name.coef = "grouptreatment", type = "estimation"))
    }
    if(fast){
        resTempo <- getInformation(out$interim.lmm, data = dataL.interim, name.coef = "grouptreatment", type = "prediction", method.prediction = "pooling2", details = TRUE)
        out$info <- c(out$info,"info.inflation" = attr(resTempo,"details")$decision$information)
        out$info <- c(out$info,"info.pooling2" = as.numeric(resTempo))
    }else{
        out$info <- c(out$info,"info.inflation" = getInformation(out$interim.lmm, data = dataL.interim, name.coef = "grouptreatment", type = "prediction", method.prediction = "inflation"))
        out$info <- c(out$info,"info.pooling2" = getInformation(out$interim.lmm, data = dataL.interim, name.coef = "grouptreatment", type = "prediction", method.prediction = "pooling2"))
    }
    out$info <- c(out$info,"rho.GS" = as.double(attr(data,"param")["rho"]), "rho" = rho)
    
    ## *** optional decision analysis (3 months, reject or not H0)
    data.decision <- data[time == decision.time]
    out$decision.ttest <- ttest(outcome ~ group, data = data.decision)
    if(fast){
        out$info <- c(out$info, "info.decision" = 1/out$decision.ttest$stderr^2)
    }else{
        out$info <- c(out$info, "info.decision" = getInformation(out$decision.ttest, data = data.decision, type = "estimation"))
    }
    
    ## *** final decision analysis (4 months, reject or not H0)
    data.final <- data[time == final.time]
    out$final.ttest <- ttest(outcome ~ group, data = data.final)
    if(fast){
        out$info <- c(out$info, "info.final" = 1/out$final.ttest$stderr^2)
    }else{
        out$info <- c(out$info, "info.final" = getInformation(out$final.ttest, data = data.final, type = "estimation"))
    }

    ## export
    return(out)
}




##----------------------------------------------------------------------
### FCT_fixedSimpleGSD.R ends here
