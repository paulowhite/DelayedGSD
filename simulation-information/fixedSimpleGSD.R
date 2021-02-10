### fixedSimpleGSD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 10 2020 (15:20) 
## Version: 
## Last-Updated: feb 10 2021 (14:25) 
##           By: Brice Ozenne
##     Update #: 267
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## path 
if(system("whoami",intern=TRUE)=="paul"){  
    setwd("~/research/SeqDesignDelayed/DelayedGSD/")
}else if(system("whoami",intern=TRUE)=="brice"){  
    setwd("~/Documents/GitHub/DelayedGSD/")
}
source("./simulation-information/FCT_fixedSimpleGSD.R")

## * R packages
library(data.table)
library(nlme)
library(pbapply)
library(ggplot2)

## * Parameters
## - 2 groups
## - 50 patients per month - max 200
## - response at 1 month
## - 1 interim analysis
## - sd = 1.2268
## - mu = c(0,0.39814)
## power.t.test(n=150,sd=1.2268,sig.level=0.05,power=0.8,type="two.sample")$delta
## power.t.test(n=100,sd=1,sig.level=0.05,power=0.8,type="two.sample")$delta

## total information (assuming sd = 1)
## 3n/(2sigma^2) = (z\alpha/2 + z\beta)^2/\delta^2
## (150)/(2*1^2)
## ((qnorm(0.975)+qnorm(0.8))/power.t.test(n=150,sd=1,sig.level=0.05,power=0.8,type="two.sample")$delta)^2

## * Simulation
cpus <- 4

cl <- snow::makeSOCKcluster(cpus)
doSNOW::registerDoSNOW(cl)
parallel::clusterExport(cl, varlist = c("analyzeData","simData"))

n.sim <- 100
ls.res <- pblapply(1:n.sim,function(iSim){    
    out <- list("0" = analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0, n.batch = 4)),
                "0.1" = analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.1, n.batch = 4)),
                "0.2" = analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.2, n.batch = 4)),
                "0.3" = analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.3, n.batch = 4)),
                "0.4" = analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.4, n.batch = 4)),
                "0.5" = analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.5, n.batch = 4)),
                "0.6" = analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.6, n.batch = 4)),
                "0.7" = analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.7, n.batch = 4)),
                "0.8" = analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.8, n.batch = 4)),
                "0.9" = analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.9, n.batch = 4)),
                "0.99" = analyzeData(simData(n = 50, sigma2 = 1.2268^2, mu0 = 0, mu1 = 0.39814, rho = 0.99, n.batch = 4))
                )
    return(as.data.table(do.call(rbind,lapply(out,"[[","info"))))
}, cl = cl)
parallel::stopCluster(cl)

## ls.res <- lapply(ls.res,"[[","info")
saveRDS(ls.res, file = paste0("simulation-information/simInfo-",n.sim,".rds"))
## ls.res <- readRDS(file = paste0("simulation-information/simInfo-",10000,".rds"))

## * Process
dt.info <- do.call(rbind,ls.res)
## setnames(dt.info, old = "rho.GS.rho", new = "rho.GS")
dt.info[,.(estimate = mean(.SD$rho), sd = sd(.SD$rho)),by="rho.GS"]
##     rho.GS    estimate          sd
##  1:   0.00 0.001433857 0.072437835
##  2:   0.10 0.097928854 0.070577009
##  3:   0.20 0.201803075 0.066222699
##  4:   0.30 0.301682046 0.062741845
##  5:   0.40 0.398217878 0.059539635
##  6:   0.50 0.499887585 0.050478667
##  7:   0.60 0.602921105 0.046484481
##  8:   0.70 0.701866189 0.034615644
##  9:   0.80 0.799848634 0.024459006
## 10:   0.90 0.900754604 0.013054795
## 11:   0.99 0.989972019 0.001297846

dtL.info <- melt(dt.info, id.vars = c("rho.GS"),
                 measure.vars = c("info.ttest","info.ttest2","info.lmm","info.inflation","info.pooling2","info.decision"),
                 variable.name = "method", value.name =  "information")
dtL.info[,method := gsub("info\\.","",method)]
dtL.info[,error := .SD$information - dtL.info[method=="decision",information], by=method]

dtL.info$method <- as.factor(dtL.info$method)
dtL.info$rho.GS <- as.factor(dtL.info$rho.GS)

## *  Display
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
## 
dtL.info[,method2 := factor(method,
                            levels = c("decision","ttest","ttest2","lmm", "inflation", "pooling2"),
                            labels = c("oracle","ttest","strategy 1","lmm", "strategy 2", "strategy 2 bis"))]

ggInfo0 <- ggplot(dtL.info[method %in% c("decision","ttest","lmm")], aes(y = information, fill = method2, x = rho.GS))
ggInfo0 <- ggInfo0 + geom_boxplot() + ylab(range(dtL.info$information))
ggInfo0 <- ggInfo0 + scale_fill_manual("", values = ggthemes::colorblind_pal()(8)[2:4]) + theme(text = element_text(size=18), legend.position="bottom")
ggInfo0 <- ggInfo0 + xlab("Correlation between X and Y") + ylab("Information")

ggError0 <- ggplot(dtL.info[method %in% c("decision","ttest2","lmm")], aes(y = error, fill = method2, x = rho.GS))
ggError0 <- ggError0 + geom_boxplot()
ggError0 <- ggError0 + scale_fill_manual("", values = ggthemes::colorblind_pal()(8)[-1]) + theme(text = element_text(size=18), legend.position="bottom")
ggError0 <- ggError0 + xlab("Correlation between X and Y") + ylab("Difference with the information at decision")

ggInfo <- ggplot(dtL.info[method %in% c("decision","ttest2","inflation")],
                 aes(y = information, fill = method2, x = rho.GS))
ggInfo <- ggInfo + geom_boxplot()
ggInfo <- ggInfo + scale_fill_manual("", values = ggthemes::colorblind_pal()(8)[-1]) + theme(text = element_text(size=18), legend.position="bottom")
ggInfo <- ggInfo + xlab("Correlation between X and Y") + ylab("Information") + coord_cartesian(ylim = c(35,65))

ggError <- ggplot(dtL.info[method %in% c("decision","ttest2","inflation")], aes(y = error, fill = method2, x = rho.GS))
ggError <- ggError + geom_boxplot()
ggError <- ggError + scale_fill_manual("", values = ggthemes::colorblind_pal()(8)[-1]) + theme(text = element_text(size=18), legend.position="bottom")
ggError <- ggError + xlab("Correlation between X and Y") + ylab("Difference with the information at decision") + coord_cartesian(ylim = c(-10,10))

ggsave(ggInfo0, filename = "./simulation-information/figures/fig-simInfo-info0.pdf", width = 10)
ggsave(ggInfo, filename = "./simulation-information/figures/fig-simInfo-info.pdf", width = 10)
ggsave(ggError0, filename = "./simulation-information/figures/fig-simInfo-error0.pdf", width = 10)
ggsave(ggError, filename = "./simulation-information/figures/fig-simInfo-error.pdf", width = 10)

## ggplot(dtL.info, aes(y = information, group = method, color = method, x = rho.GS)) + geom_smooth()

dtL.info[rho.GS==0,.(rep = .N, bias = mean(error), sd = sd(error), mse = mean(error^2)), by = "method"]
##       method rep         bias       sd        mse
## 1:     ttest 100 -16.78081788 2.466855 287.620368
## 2:    ttest2 100  -0.09778406 3.068602   9.331718
## 3:       lmm 100 -16.71651461 2.491968 285.589667
## 4: inflation 100  -0.09646969 3.072653   9.356089
## 5:  pooling2 100  -0.09396540 3.073973   9.363647
## 6:  decision 100   0.00000000 0.000000   0.000000

dtL.info[rho.GS==0.5,.(rep = .N, bias = mean(error), sd = sd(error), mse = mean(error^2)), by = "method"]
##       method rep        bias       sd        mse
## 1:     ttest 100 -17.0578048 2.198750 295.754859
## 2:    ttest2 100  -0.2741210 2.443929   5.988202
## 3:       lmm 100 -14.0748119 2.437889 203.984201
## 4: inflation 100  -0.3450121 2.489230   6.253336
## 5:  pooling2 100  -0.3640734 2.554245   6.591474
## 6:  decision 100   0.0000000 0.000000   0.000000

dtL.info[rho.GS==0.7,.(rep = .N, bias = mean(error), sd = sd(error), mse = mean(error^2)), by = "method"]
##       method rep         bias       sd        mse
## 1:     ttest 100 -16.73745706 2.449738 286.083672
## 2:    ttest2 100   0.03366496 3.227454  10.313430
## 3:       lmm 100 -10.09436625 2.290147 107.088554
## 4: inflation 100   0.10966948 2.794487   7.743096
## 5:  pooling2 100   0.14585720 2.768736   7.610516
## 6:  decision 100   0.00000000 0.000000   0.000000

dtL.info[rho.GS==0.99,.(rep = .N, bias = mean(error), sd = sd(error), mse = mean(error^2)), by = "method"]
##       method rep         bias        sd         mse
## 1:     ttest 100 -16.40320856 2.6874014 276.2151560
## 2:    ttest2 100   0.71084670 3.3091392  11.3462009
## 3:       lmm 100  -0.52974115 0.7216799   0.7962393
## 4: inflation 100  -0.02391864 0.7250516   0.5210149
## 5:  pooling2 100  -0.02981053 0.7226524   0.5178929
## 6:  decision 100   0.00000000 0.0000000   0.0000000



######################################################################
### fixedSimpleGSD.R ends here
