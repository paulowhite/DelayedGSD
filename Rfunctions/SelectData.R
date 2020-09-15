### SelectData.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 28 2020 (16:33) 
## Version: 
## Last-Updated: Aug 28 2020 (17:00) 
##           By: Paul Blanche
##     Update #: 28
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
    N.fw <- length(grep("t",names(d)))
    dd <- d[d$t1 <= t - Delta.t,]
    for(j in 2:N.fw){
        whereNA <- which(dd[,paste0("t",j)] > t - Delta.t)
        dd[whereNA,paste0("t",j)] <- NA
        dd[whereNA,paste0("X",j)] <- NA
    }
    dd
}


#----------------------------------------------------------------------
### SelectData.R ends here
