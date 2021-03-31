### SelectData.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 28 2020 (16:33) 
## Version: 
## Last-Updated: mar 26 2021 (15:04) 
##           By: Brice Ozenne
##     Update #: 43
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


##' @title Select available longitudnal data for analysis at a specific follow-up time.
##' @description yy
##' @param d data set generated with GenData function.
##' @param Delta.t
##' @param t time at which we want the available data.
##' @details xxx
##' @return input dataset  with less rows and NA where appropriate.
##' @author Paul Blanche
##' 
##' @examples
##' x <- GenData(n=50)
##' head(x$d,n=20)
##' tail(x$d)
##' y <- SelectData(x$d,t=4)
##' head(y)
##' tail(y)
##' PlotProgress(x$d,at=4)
##'
##' #----- select when half of the subjects have one follow-up measuement---
##' # which depends on accrual rate (ar) and time to process data (Delta.t)
##'
##' x <- GenData(n=35)
##' theDelta.t <- 0.500001
##' thear <- 10
##' thet <- x$d$t2[ceiling(nrow(x$d)/2) + ceiling(thear*theDelta.t)]
##' PlotProgress(x$d,at=thet)
##' y <- SelectData(x$d,t=thet)
##' y
##'
##' 
SelectData <- function(d,t,Delta.t=0.500001){

    ## find columns in the dataset d corresponding to the outcome and the time
    col.outcome <- grep("X",names(d), value = TRUE)
    col.times <- grep("t",names(d), value = TRUE)
    n.times <- length(col.times)

    ## add indicator of missing values
    col.missing <- gsub("^t","missing",col.times)
    test.na <- 1.0 * is.na(d[,col.outcome,drop=FALSE]) ## 1.0 * to convert to numeric
    colnames(test.na) <- col.missing
    d.augm <- cbind(d, test.na)

    ## subset according to time
    dd <- d.augm[which(d.augm$t1 <= t - Delta.t),,drop=FALSE]

    ## add missing value according to time as the baseline can be observed but some of the follow-up value may be not
    for(j in 2:n.times){
        whereNA <- which(dd[,paste0("t",j)] > t - Delta.t)
        dd[whereNA,col.missing[[j]]] <- -1
        dd[whereNA,paste0("t",j)] <- NA
        dd[whereNA,paste0("X",j)] <- NA
    }

    ## keep arguments
    attr(dd,"t") <- t
    attr(dd,"Delta.t") <- Delta.t
    return(dd)
}


#----------------------------------------------------------------------
### SelectData.R ends here
