library(DelayedGSD)
## * read all results
## library(data.table)
## res2stage <- readRDS(file.path("Simulations","brice-biostat-cluster","Results-built","res2stage.rds"))
## res2stage[missing==TRUE&binding==TRUE&fixC==TRUE&ar==5&hypo=="power"&method==1&seed==65836753]
## * example 1

method <- 1
binding <- TRUE
fixC <- TRUE
seed <- 65836753

e.plannedB <- CalcBoundaries(kMax = 2,  
                             alpha = 0.025, 
                             beta = 0.2,  
                             InfoR.i = c(0.482, 1.00),  
                             InfoR.d = c(0.52,1),  
                             rho_alpha = 2,  
                             rho_beta = 2,  
                             method = method,  
                             cNotBelowFixedc = fixC,
                             bindingFutility= binding,
                             delta = 0.6)

df.sim <- GenData(n = 486, 
                  N.fw = 2,
                  rand.block = c(1,1,0,0),
                  allsd = c(2.5,2.1,2.4),
                  mean0 = c(10,0,0),
                  delta = c(0,0.3,0.6),
                  ar = (0.86*2)*2*5,
                  cor.01.1 = -0.15,
                  cor.ij.1 = 0.68,
                  cor.0j.1 = -0.27,
                  seed = seed,
                  MissProb = matrix(c(0, 0, 0, 0), 
                                    nrow = 2, 
                                    ncol = 2, 
                                    dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing")) 
                                    ) ,
                  DigitsOutcome = 2,
                  TimeFactor = 14,
                  DigitsTime = 0)$d
thets <- c(218,218,218)
nGSD <- c(486,486,486)

## ** interim
## SelectData(df.sim,t=thets[method])
df.simI <- SelectData(df.sim,t=thets[method])
e.lmmI <- analyzeData(df.simI,
                      ddf = "nlme", data.decision = sum(df.sim$t1 <= thets[method] + 1.50001*14), getinfo = TRUE, trace = TRUE)
e.GSDI <- update(e.plannedB, delta = e.lmmI, trace = TRUE)

e.GSDI$ck
e.GSDI$ck.unrestricted
## [1] 1.485608

## ** decision
e.lmmD <- analyzeData(df.sim[which(df.sim$t1 <= thets[method] + 1.50001*14),],
                      ddf = "nlme", getinfo = TRUE, trace = TRUE)
if(e.GSDI$conclusion["interim",1]=="stop"){
    e.GSDD <- update(e.GSDI, delta = e.lmmD, trace = FALSE)
}else{
    e.GSDD <- update(e.GSDI, delta = e.lmmD, k = 1, type.k = "decision", trace = FALSE)
}
e.GSDD$ck
e.GSDD$ck.unrestricted
## [1] 1.49773

## ** final
if(e.GSDD$conclusion["interim",1]=="continue"){
    e.lmmF <- analyzeData(df.sim[1:nGSD[method],],
                          ddf = "nlme", getinfo = TRUE, trace = TRUE)
    e.GSDF <- update(e.GSDD, delta = e.lmmF, trace = FALSE)
    summary(e.GSDF)
    ## confint(debug.GSD)

    seqC <- c(e.GSDF$ck.unrestricted, e.GSDF$uk[2])
    seqDelta <- unique(round(sort(c(seq(0,1.2, length.out = 10), confint(e.GSDF)[2,"lower"], confint(e.GSDF)[2,"upper"])),2))
    ## seqZ <- sort(c(qnorm(0.975),seq(-5,5, length.out = 4)))
    seqZ <- unique(sort(c(seqC, confint(e.GSDF)[1,"statistic"], qnorm(0.975),seq(0,5, length.out = 10),
                          seq(seqC[2],qnorm(0.975),length.out=4))))


    grid <- expand.grid(stage = 1:2,                        
                        z = seqZ,
                        delta = seqDelta)
    grid$reject <- (grid$stage==1)*(grid$z>=max(seqC[1], qnorm(0.975))) + (grid$stage==2)*(grid$z>=max(seqC[2], qnorm(0.975)))
    grid$p.value <- as.numeric(NA)
    NROW(grid)

    grid0 <- grid[grid$reject==0,]
    grid1 <- grid[grid$reject==1,]
    grid <- rbind(grid0[order(grid0$stage,grid0$z),],
                  grid1[order(-grid1$stage,grid1$z),])
    n.grid <- NROW(grid)

    for(iGrid in 1:n.grid){ ## iGrid <- 386 375
        cat("*")
        ## which
        iStage <- grid[iGrid,"stage"] ## iStage <- 1
        iZ <- grid[iGrid,"z"] ## iZ <- seqZ[7]
        iDelta <- grid[iGrid,"delta"] ## iDelta <- 0

        if(iStage == 1){
            if(iZ >= seqC[1] && iZ < qnorm(0.975)){
                iZ <- seqC[1]
            }
            iEstimate <- iZ / sqrt(12.58797814)
            grid[iGrid,"p.value"] <- FinalPvalue(Info.d = c(12.58797814),  
                                                 Info.i = c(10.60271846),
                                                 ck = c(1.49772992),
                                                 lk = c(0.24335814),  
                                                 uk = c(2.5458844),  
                                                 kMax = 2, 
                                                 delta = iDelta,  
                                                 estimate = iEstimate,
                                                 method = 1,
                                                 bindingFutility = TRUE,
                                                 cNotBelowFixedc = TRUE)
        }else if(iStage == 2){
            if(iZ >= seqC[2] && iZ < qnorm(0.975)){
                iZ <- seqC[2]
            }
            iEstimate <- iZ / sqrt(20.74866926) ## 1.9927 / sqrt(20.74866926)
            grid[iGrid,"p.value"] <- FinalPvalue(Info.d = c(12.58797814),  
                                                 Info.i = c(10.60271846, 20.74866926),
                                                 ck = c(1.49772992),
                                                 lk = c(0.24335814, 1.9961649),  
                                                 uk = c(2.5458844, 1.9961649),  
                                                 kMax = 2, 
                                                 delta = iDelta,  
                                                 estimate = iEstimate, ## 0.4374693
                                                 method = 1,
                                                 bindingFutility = TRUE,
                                                 cNotBelowFixedc = TRUE)
        }
        
    }

    UstageZ <- unique(grid[,c("reject","stage","z")])
    UstageZ$index <- 1:NROW(UstageZ)

    gridI <- merge(UstageZ,grid, by = c("reject","stage","z"))
    gridI <- gridI[order(gridI$index),]

    ## plot 1
    for(iD in 1:length(seqDelta)){ ## iD <- 1
        iDelta <- seqDelta[iD]
        iGrid <- gridI[gridI$delta == iDelta,]
        if(iD==1){
            plot(x = iGrid$index, y = iGrid$p.value, col = iD, pch = iGrid$stage, axes = FALSE, type = "b",
                 xlab = "z", ylab = "Pr(Z>z|theta)")
        }
        ## else{
            lines(x = iGrid$index, y = iGrid$p.value, col = iD, pch = iGrid$stage, type = "b")
        ## }
        axis(1, at = 1:NROW(UstageZ), labels = round(UstageZ$z,2))
        axis(2, las = 2)
    }
    ## abline(v = which(UstageZ$stage==2)[1]-0.5, lty = 2)
    nZ.D1 <- sum(seqZ<max(seqC[1],qnorm(0.975)))
    nZ.F2 <- nZ.D1 + length(seqZ)
    nZ.D2 <- NROW(UstageZ)

    index.obs <- which((UstageZ$stage==2) & (UstageZ$z==confint(e.GSDF)[1,"statistic"]))
    abline(h = c(0.025,0.975), lty = 2, col = "gray")
    abline(v = c(nZ.D1,nZ.F2,nZ.D2)+0.5, lty = 2)
    abline(v = index.obs, lwd = 2, lty = 2)
    
    legend("bottomleft", legend = paste0("theta=",seqDelta), col = 1:length(seqDelta), pch = 20)

    ## confint(e.GSDF)[1,"statistic"]
    ## To do:
    ##   - print CI values in pot
    ##   - shaded area stages 1,2

    ## plot 2
    grid.subset <- gridI[gridI$index == index.obs,]
    plot(x = grid.subset$delta, y = grid.subset$p.value, type = "b",
         xlab = "theta", ylab = "Pr(Z>z|theta)")
    abline(h = c(0.025,0.975), lty = 2, col = "gray")
    abline(v = c(confint(e.GSDF)[2,"lower"],confint(e.GSDF)[2,"upper"]), lty = 2)
    
    
}


