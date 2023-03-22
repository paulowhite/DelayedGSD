### figure-pvalue-correction.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 10 2023 (13:16) 
## Version: 
## Last-Updated: mar 21 2023 (18:44) 
##           By: Brice Ozenne
##     Update #: 26
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(DelayedGSD)
library(RColorBrewer)
export <- FALSE

## * 2 stages

## ** graphical display - no fix C
if(export){
    pdf("figures/illustration-pvalue-2stage.pdf", width = 5)
}
res2stage <- gridFinalPvalue(list(Info.d = c(12.58797814),
                                  Info.i = c(10.60271846, 20.74866926),
                                  ck = 1.49772992, ck.unrestricted = 1.49772992,
                                  lk = c(0.24335814, 1.9961649),
                                  uk = c(2.5458844, 1.9961649),
                                  kMax = 2,
                                  reason.interim = c("no boundary crossed",NA),
                                  method = 1,
                                  bindingFutility = TRUE, 
                                  cNotBelowFixedc = FALSE),
                             continuity.correction = 0,
                             xlim = c(0.7,2.2))
if(export){
    dev.off()
}

## ** graphical display - fix C
res2stage.fixC <- vector(mode = "list", length = 3)

if(export){
    pdf("figures/illustration-pvalue-2stage-fixC.pdf", width = 12)
}
par(mfrow = c(1,3), mar = rep(4,4))
for(iC in 0:2){ ## iC <- 0
    iTitle <- switch(as.character(iC),
                     "0" = "no correction",
                     "1" = "continuity correction \n (add P)",
                     "2" = "continuity correction \n (shift stat)",
                     )
    set.seed(10)
    res2stage.fixC[[iC+1]] <- gridFinalPvalue(list(Info.d = c(12.58797814),
                                                   Info.i = c(10.60271846, 20.74866926),
                                                   ck = 1.959964, ck.unrestricted = 1.49772992,
                                                   lk = c(0.24335814, 1.9961649),
                                                   uk = c(2.5458844, 1.9961649),
                                                   kMax = 2,
                                                   reason.interim = c("no boundary crossed",NA),
                                                   method = 1,
                                                   bindingFutility = TRUE, 
                                                   cNotBelowFixedc = TRUE),
                                              continuity.correction = iC,
                                              xlim = c(0.7,2.2), title = iTitle)
}
if(export){
    dev.off()
}

## * 3 stages

## ** graphical display - no fix C
bound3s <- CalcBoundaries(kMax = 3,  
                          alpha = 0.025, 
                          beta = 0.2,  
                          InfoR.i = c(0.40, 0.7, 1.00),  
                          InfoR.d = c(0.45, 0.75, 1),  
                          rho_alpha = 2,  
                          rho_beta = 2,  
                          method = 1,  
                          cNotBelowFixedc = FALSE,
                          bindingFutility= TRUE,
                          delta = 0.6)

## Planned boundaries: 
## stage  F-bound E-bound C-bound alpha-spent beta-spent
##     1 -0.01468 2.65207 1.39026     0.00400      0.032
##     2  1.07723 2.32431 1.74696     0.01225      0.098
##     3                  2.03284     0.02500      0.200

## Planned information: 
## stage  Interim (%) Decision  (%)
##     1  9.37888 0.4 10.55124 0.45
##     2 16.41305 0.7 17.58541 0.75
##     3              23.44721 1.00

if(export){
    pdf("figures/illustration-pvalue-3stage.pdf", width = 5)
}
res3stage <- gridFinalPvalue(bound3s,
                             continuity.correction = 0,
                             xlim = c(0.7,3.3))
if(export){
    dev.off()
}


## ** graphical display - no fix C
bound3s.fixC <- CalcBoundaries(kMax = 3,  
                               alpha = 0.025, 
                               beta = 0.2,  
                               InfoR.i = c(0.40, 0.7, 1.00),  
                               InfoR.d = c(0.45, 0.75, 1),  
                               rho_alpha = 2,  
                               rho_beta = 2,  
                               method = 1,  
                               cNotBelowFixedc = TRUE,
                               bindingFutility= TRUE,
                               delta = 0.6)

res3stage.fixC <- vector(mode = "list", length = 3)
## Planned boundaries: 
## stage  F-bound E-bound C-bound alpha-spent beta-spent
##     1 -0.01468 2.65207 1.95996     0.00400      0.032
##     2  1.07723 2.32431 1.95996     0.01225      0.098
##     3                  2.03284     0.02500      0.200

## Planned information: 
## stage  Interim (%) Decision  (%)
##     1  9.37888 0.4 10.55124 0.45
##     2 16.41305 0.7 17.58541 0.75
##     3              23.44721 1.00

if(export){
    pdf("figures/illustration-pvalue-3stage-fixC.pdf", width = 12)
}
par(mfrow = c(1,3), mar = rep(4,4))
for(iC in 0:2){ ## iC <- 0
    iTitle <- switch(as.character(iC),
                     "0" = "no correction",
                     "1" = "continuity correction \n (add P)",
                     "2" = "continuity correction \n (shift stat)",
                     )
    set.seed(10)
    res3stage.fixC[[iC+1]] <- gridFinalPvalue(bound3s.fixC,
                                              continuity.correction = iC,
                                              xlim = c(0.65,3.3),
                                              title = iTitle,
                                              digits = 3)
}
if(export){
    dev.new()
}

## ** calculation of the p-value
## CalcBoundaries(kMax = 3,  
##                alpha = 0.025, 
##                beta = 0.2,  
##                InfoR.i = c(0.40, 0.7, 1.00),  
##                InfoR.d = c(0.45, 0.75, 1),  
##                rho_alpha = 2,  
##                rho_beta = 2,  
##                method = 1,  
##                cNotBelowFixedc = FALSE,
##                bindingFutility= TRUE,
##                delta = 0.6)
## alpha-spent
##     0.00400
##     0.01225
##     0.02500

calcP_3stage <- function(z, k, continuity.correction){
    if(k==1){
        out <- FinalPvalue(Info.d = c(10.55124),
                           Info.i = c(9.37888),
                           ck = c(1.959964), ck.unrestricted = c(1.39026),
                           lk = c(-0.01468),
                           uk = c(2.65207),
                           kMax = 3,
                           delta = 0, 
                           estimate = z / sqrt(10.55124),
                           reason.interim = c("efficacy",NA,NA),
                           method = 1,
                           bindingFutility = TRUE, 
                           cNotBelowFixedc = TRUE,
                           continuity.correction = continuity.correction)
    }else if(k==2){
        out <- FinalPvalue(Info.d = c(10.55124, 17.58541),
                           Info.i = c(9.37888, 16.41305),
                           ck = c(1.959964,1.959964), ck.unrestricted = c(1.39026,1.74696),
                           lk = c(-0.01468, 1.07723),
                           uk = c(2.65207, 2.32431),
                           kMax = 3,
                           delta = 0, 
                           estimate = z / sqrt(17.58541),
                           reason.interim = c("no boundary crossed","efficacy",NA),
                           method = 1,
                           bindingFutility = TRUE, 
                           cNotBelowFixedc = TRUE,
                           continuity.correction = continuity.correction)
    }else if(k==3){
        out <- FinalPvalue(Info.d = c(10.55124, 17.58541),
                           Info.i = c(9.37888, 16.41305, 23.44721),
                           ck = c(1.959964,1.959964), ck.unrestricted = c(1.39026,1.74696),
                           lk = c(-0.01468, 1.07723, 2.03284),
                           uk = c(2.65207, 2.32431, 2.03284),
                           kMax = 3,
                           delta = 0, 
                           estimate = z / sqrt(23.44721),
                           reason.interim = c("no boundary crossed","no boundary crossed",NA),
                           method = 1,
                           bindingFutility = TRUE, 
                           cNotBelowFixedc = TRUE,
                           continuity.correction = continuity.correction)
    }

    return(data.frame(p.value = out, z = z, k = k, continuity.correction = continuity.correction))
}

df.p3s <- NULL
for(iC in 0:2){ ## iC <- 0
    df.p3s <- rbind(df.p3s,
                  fut.1_1 = calcP_3stage(z = -3, k = 1, continuity.correction = iC),
                  fut.1_2 = calcP_3stage(z = 0, k = 1, continuity.correction = iC),
                  fut.1_Bl = calcP_3stage(z = 1.39026, k = 1, continuity.correction = iC),
                  fut.1_Bu = calcP_3stage(z = 1.95, k = 1, continuity.correction = iC),
                  fut.2_1 = calcP_3stage(z = -3, k = 2, continuity.correction = iC),
                  fut.2_2 = calcP_3stage(z = 1, k = 2, continuity.correction = iC),
                  eff.2_Bl = calcP_3stage(z = 1.74696, k = 2, continuity.correction = iC),
                  eff.2_Bu = calcP_3stage(z = 1.95, k = 2, continuity.correction = iC),
                  fut.3_1 = calcP_3stage(z = -3, k = 3, continuity.correction = iC),
                  fut.3_2 = calcP_3stage(z = 1, k = 3, continuity.correction = iC),
                  eff.3_B = calcP_3stage(z = 2.03284, k = 3, continuity.correction = iC),
                  eff.3_1 = calcP_3stage(z = 3, k = 3, continuity.correction = iC),
                  eff.3_2 = calcP_3stage(z = 6, k = 3, continuity.correction = iC),
                  eff.2_B = calcP_3stage(z = 1.959964, k = 2, continuity.correction = iC),
                  eff.2_1 = calcP_3stage(z = 2.5, k = 2, continuity.correction = iC),
                  eff.2_1 = calcP_3stage(z = 3.25, k = 2, continuity.correction = iC),
                  eff.2_1 = calcP_3stage(z = 4, k = 2, continuity.correction = iC),
                  eff.2_2 = calcP_3stage(z = 5, k = 2, continuity.correction = iC),
                  eff.1_B = calcP_3stage(z = 1.959964, k = 1, continuity.correction = iC),
                  eff.1_1 = calcP_3stage(z = 2.5, k = 1, continuity.correction = iC),
                  eff.1_1 = calcP_3stage(z = 3.25, k = 1, continuity.correction = iC),
                  eff.1_1 = calcP_3stage(z = 4, k = 1, continuity.correction = iC),
                  eff.1_2 = calcP_3stage(z = 5, k = 1, continuity.correction = iC)
                  )
}
rownames(df.p3s) <- NULL
df.p3s$p.display <- paste0(round(100*df.p3s$p.value,2),"%")
## df.p[df.p$k==1 & df.p$z == 1.95,"z.display"] <- 1.9


## ** graphical display
vec.col <- tail(brewer.pal(n = 8, name =  "PuBu"),5)


pdf("figures/illustration3s-pvalue-correction.pdf", width = 12)
par(mfrow = c(1,3), mar = rep(4,4))
displayPstage(data = df.p3s[df.p3s$continuity.correction==0,],
              xlim = c(0.7,3.3), title = "no correction", col = vec.col)
displayPstage(data = df.p3s[df.p3s$continuity.correction==1,],
              xlim = c(0.7,3.3), title = "continuity correction \n (add P)", col = vec.col)
displayPstage(data = df.p3s[df.p3s$continuity.correction==2,],
              xlim = c(0.7,3.3), title = "continuity correction \n (shift stat)", col = vec.col)
dev.off()





##----------------------------------------------------------------------
### figure-pvalue-correction.R ends here
