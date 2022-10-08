### FCT.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  7 2022 (14:40) 
## Version: 
## Last-Updated: okt  7 2022 (15:41) 
##           By: Brice Ozenne
##     Update #: 42
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

exportGSD <- function(object,
                      export.statistic = TRUE,
                      export.ML = TRUE,
                      export.MUE = TRUE,
                      export.info = TRUE,
                      export.predinfo = TRUE,
                      export.boundary = TRUE,
                      export.decision = TRUE){


    ## ** initialize
    out <- data.frame(statistic = NA,
                      estimate_ML = NA,
                      se_ML = NA,
                      p.value_ML = NA,
                      lower_ML = NA,
                      upper_ML = NA,
                      estimate_MUE = NA,
                      p.value_MUE = NA,
                      lower_MUE = NA,
                      upper_MUE = NA,
                      info = NA,
                      infoPC = NA,
                      info.pred = NA,
                      infoPC.pred = NA,
                      uk = NA,
                      lk = NA,
                      ck = NA,
                      decision = NA,
                      reason = NA)
    
    ## ** check user input
    if(identical(object,NA)){
       return(out) 
    }else if(!inherits(object,"delayedGSD")){
        stop("Argument \'object\' inherits from \"delayedGSD\". \n")
    }

    ## ** extract information
    stage <- object$stage
    if(stage$type %in% "interim"){
        object.confint <- confint(object)
    }else  if(stage$type %in% c("decision","final")){
        object.confint <- confint(object, method = c("ML","MUE"))
    }
    object.confint <- object.confint[object.confint$stage==stage$k,,drop=FALSE]
    
    object.info <- coef(object, type = "information")
    object.info <- object.info[object.info$stage==stage$k,,drop=FALSE]

    object.boundary <- coef(object, type = "boundary")
    object.boundary <- object.boundary[object.boundary$stage==stage$k,,drop=FALSE]

    object.decision <- coef(object, type = "decision")

    ## ** fill output
    if(export.statistic){
        out$statistic <- object.confint[1,"statistic"]
    }
    if(export.ML){
        if(stage$type %in% "interim"){
            out$estimate_ML <- object.confint[,"estimate"]
            out$se_ML <- object.confint[,"se"]
        }else if(stage$type %in% c("decision","final")){
            out$estimate_ML <- object.confint[object.confint$method == "ML","estimate"]
            out$se_ML <- object.confint[object.confint$method == "ML","estimate"]
            out$p.value_ML <- object.confint[object.confint$method == "ML","p.value"]
            out$lower_ML <- object.confint[object.confint$method == "ML","lower"]
            out$upper_ML <- object.confint[object.confint$method == "ML","upper"]
        }
    }

    if(export.MUE){
        if(stage$type %in% c("decision","final")){
            out$estimate_MUE <- object.confint[object.confint$method == "MUE","estimate"]
            out$p.value_MUE <- object.confint[object.confint$method == "MUE","p.value"]
            out$lower_MUE <- object.confint[object.confint$method == "MUE","lower"]
            out$upper_MUE <- object.confint[object.confint$method == "MUE","upper"]
        }
    }

    if(export.info){
        if(stage$type %in% c("interim","final")){
            out$info <- object.info[,"Interim"]
            out$infoPC <- object.info[,"Interim.pc"]
        }else if(stage$type == "decision"){
            out$info <- object.info[,"Decision"]
            out$infoPC <- object.info[,"Decision.pc"]
        }
    }

    if(export.predinfo){
        if(stage$type %in% "interim"){
            out$info.pred <- object.info[,"Decision"]
            out$infoPC.pred <- object.info[,"Decision.pc"]
        }
    }

    if(export.boundary){
        if(stage$type %in% "interim"){
            out$uk <- object.boundary[,"Ebound"]
            out$lk <- object.boundary[,"Fbound"]
        }else if(stage$type %in% c("decision","final")){
            out$ck <- object.boundary[,"Cbound"]
        }
    }

    if(export.decision){
        out$decision <- object.decision["decision",NCOL(object.decision)]
        out$reason <- object.decision["comment",NCOL(object.decision)]
    }


    ## ** export
    return(out)

}


##----------------------------------------------------------------------
### FCT.R ends here
