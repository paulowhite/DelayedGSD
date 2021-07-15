## * plot.delayedGSD (documentation)
#' @title Display boundaries over information
#'
#' @param x object from CalcBoundaries
#' @param y not used. For compatibility with the generic method.
#' @param type type of boundaries to plot (Z-statistic (Z), p-value (P), effect (E))
#' @param Itype information scale on x-axis (rate or absolute (abs))
#' @param main to specify a plot title
#' @param xlim range of the x axis
#' @param ylim range of the y axis
#' @param legend should a caption be added?
#' @param legend.ncol number of column in the caption
#' @param legend.cex character expansion factor for the text in the legend.
#' @param add.arrows [logical] should arrows indicating the ordering of the space be added?
#' @param sep.arrows [numeric vector of length 2] space between the critical point and the arrow (below and above)
#' @param size.arrows [numeric, >0] width of the arrow.
#' @param lwd.arrows [numeric, >0] thickness of the arrow.
#' @param ... not used. For compatibility with the generic method.
#'
#' @examples
#' bnds <- CalcBoundaries(kMax=3, InfoR.i = c(0.5,0.75,1), InfoR.d = c(0.55,0.8))
#' plot(bnds, type = "Z")
#' plot(bnds, type = "Z", add.arrows = TRUE)
#' plot(bnds, type = "Z", add.arrows = TRUE, lwd.arrow=2,ylim=c(0.35,2.6), xlim=c(0.5,1.1),legend=FALSE)
#' plot(bnds, type = "E")
#' plot(bnds, type = "P")

## * plot.delayedGSD (code)
#' @export
plot.delayedGSD <- function(x,
                            y,
                            type="Z",   
                            Itype="rate",
                            planned=TRUE,
                            main = NULL,   
                            xlim = NULL, 
                            ylim = NULL, 
                            legend=TRUE,  
                            legend.ncol=1,
                            legend.cex=1,
                            add.arrows=FALSE,
                            sep.arrow=c(0.2,0.1),
                            size.arrow=0.1,
                            lwd.arrow=1,
                            ...){

                                        # {{{ Check user input
    if(!missing(y)){
        stop("Arugment \'y\' not used. \n")
    }
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    type <- match.arg(type, c("Z","E","P"))
    Itype <- match.arg(Itype, c("rate","abs"))
                                        # }}}
    
                                        # {{{ Preliminaries

    specdec <- function(x,k){ format(round(x,k),nsmall=k)}

    ## extract from object
    alpha <- x$alpha
    kMax <- x$kMax
    k <- x$stage$k
    type.k <- x$stage$type
    test.planning <- identical(planned,"only") || (type.k=="planning")

    ls.info <- getInformation(x, planned = planned)
    Info.i <- ls.info$Info.i
    Info.d <- ls.info$Info.d
    Info.max <- ls.info$Info.max
    uk <- ls.info$uk
    lk <- ls.info$lk
    ck <- ls.info$ck
    if(!test.planning){
        delta <- switch(type,
                        "Z" = ls.info$delta[,"statistic"],
                        "E" = ls.info$delta[,"estimate"],
                        "P" = ls.info$delta[,"p.value"])
    }else{
        delta <- NULL
    }
    
    if(is.null(main)){
        if(test.planning){
            main <- "Planning"
        }else if(type.k=="final"){
            main <- paste0("Final analysis (stage ",k,")")
        }else if(type.k=="decision"){
            main <- paste0("Decision analysis at stage ",k)
        }else if(type.k=="interim"){
            main <- paste0("Interim analysis at stage ",k)
        }
    }

    ## prepare
    if(Itype=="rate"){
        Ival <- Info.i/Info.max
        Idval <- Info.d/Info.max
        xlb <- "Information fraction"
    } else if(Itype=="abs"){
        Ival <- Info.i
        Idval <- Info.d
        xlb <- "Information"
    } 
    
    if(type=="P"){
        mydigits <- 3
    }else{
        mydigits <- 2
    }
                                        # }}}
                                        # {{{ Values to plot
                                        # default (if type="Z")
    xu <- Ival
    xd <- Idval
    yu <- uk
    yl <- lk
    yc <- ck
    yh <- qnorm(1-alpha)
    whereleg <- "bottomleft"
    ylab <- "Stopping boundary (Z-statistic)"
 
    if(type=="P"){
        whereleg <- "topleft"
        ylab <- "Stopping boundary (P-value)"
        yu <- 1-pnorm(uk)
        yl <- 1-pnorm(lk)
        yc <- 1-pnorm(ck)
        yh <- alpha
    }
    if(type=="E"){
        ylab <- "Stopping boundary (Effect estimate)"
        yu <- uk/sqrt(Info.i)
        yl <- lk/sqrt(Info.i)
        yc <- ck/sqrt(Info.d)
        yh <- qnorm(1-alpha)/sqrt(Info.i[length(Info.i)])
    }
    if(is.null(ylim)){
        ylim <- switch(type,
                       "P" = c(0,1),
                       "E" = c(min(c(yl,delta))-0.5,max(c(yu,delta))+0.5),
                       "Z" = c(min(c(-1,yl,delta)),max(c(3,yl,delta)))
                       )
    }
    if(is.null(xlim)){
        xlim <- c(0,max(xu))
    }
                                        # }}}
                                        # {{{ Plot
    plot(xu,yu,type="l",lty=2,lwd=2,ylim=ylim,col="green3",xlab=xlb,ylab=ylab,axes=FALSE,
         xlim=xlim,main=main)
    lines(xu,yl,col="red",lwd=2,lty=2)
    points(xd,yc,col="black",pch=19,cex=1.5)
    points(xu,yl,col="red",pch=21,bg="red",cex=1.2)
    points(xu,yu,col="green3",pch=21,bg="green3",cex=1.2)
    if(!is.null(delta)){
        xdelta <- xu[1:k]
        if(type.k=="decision"){
            xdelta <- c(xdelta,xd[k])
        }
        ## lines(c(0,xdelta),c(0,delta),col="purple",lwd=2,lty=3)
        points(xdelta,delta,col="purple",pch=22,bg="purple",cex=1.2)
        text(x=xdelta,y=delta,labels=specdec(delta,k=mydigits),col="purple",pos=1)
        
    }
                                        #---
    text(x=xd,y=yc,labels=specdec(yc,k=mydigits),col="black",pos=1)
    text(x=xu,y=yl,labels=specdec(yl,k=mydigits),col="red",pos=1)
    text(x=xu,y=yu,labels=specdec(yu,k=mydigits),col="green3",pos=3)    
                                        #---   
    abline(h=yh,lty=2,col="grey")
    text(x=0,y=yh,labels=specdec(yh,k=mydigits),col="grey",pos=1)        
                                        # axes
    axis(1,at=c(0,xu,xd),labels=format(c(0,xu,xd),digits=2))
    axis(2,at=c(yl,yu,ylim),
         labels=format(c(yl,yu,ylim),digits=2),las=2)
                                        # color areas
    colvectgreen <- col2rgb("green3")/255
    colvectred <- col2rgb("red")/255
    myrgbcolgreen <- rgb(red=colvectgreen[1],
                         green=colvectgreen[2],
                         blue=colvectgreen[3],
                         alpha=0.25)
    myrgbcolred <- rgb(red=colvectred[1],
                       green=colvectred[2],
                       blue=colvectred[3],
                       alpha=0.25)
    polygon(x=c(0,xu,rev(xu),0),
            y=c(yl[1],yl,rep(ifelse(type=="P",max(ylim),min(ylim)),length(yu)+1)),
            col=myrgbcolred,border=NA)
    polygon(x=c(0,xu,rev(xu),0),
            y=c(yu[1],yu,rep(ifelse(type!="P",max(ylim),min(ylim)),length(yu)+1)),
            col=myrgbcolgreen,border=NA)
                                        # legend   
    if(legend){
        if(length(legend.cex)==1){
            legend.cex <- rep(legend.cex,2)
        }
        
        text.legend <- c("Stopping bound for efficacy",
                         "Stopping bound for futility",
                         "Critical value fixed design",
                         "Decision boundary")
        lty.legend  <- c(2,2,2,NA)
        pch.legend  <- c(21,21,NA,19)
        col.legend  <- c("green3","red","grey","black")
        pt.bg.legend <- c("green3","red",NA,NA)
        if(!is.null(delta)){
            text.legend <- c(text.legend,"Current estimate")
            lty.legend <- c(lty.legend,NA)
            pch.legend <- c(pch.legend,22)
            col.legend <- c(col.legend,"purple")
            pt.bg.legend <- c(pt.bg.legend,"purple")
        }
        
        legend(whereleg,
               legend = text.legend,
               lty = lty.legend,
               pch = pch.legend,
               col = col.legend,
               pt.bg = pt.bg.legend,
               bg="white",
               cex=legend.cex[1],
               pt.cex=legend.cex[2], 
               ncol = legend.ncol)
    }
    if(add.arrows){
        n.decision <- length(ck)
        for(iDecision in 1:n.decision){
            arrows(x0=Idval[iDecision],y0=ylim[1],x1=Idval[iDecision],y1=yc[iDecision]-sep.arrow[1],lty=1,length=size.arrow,lwd=lwd.arrow)
            arrows(x0=Idval[iDecision],y0=yc[iDecision]+sep.arrow[2],x1=Idval[iDecision],y1=ylim[2],lty=1,length=size.arrow,lwd=lwd.arrow)
            if(iDecision+1<=n.decision){
                arrows(x0=Idval[iDecision],y0=yc[iDecision]-sep.arrow[1],x1=Idval[iDecision+1],y1=ylim[1],lty=1,length=size.arrow,lwd=lwd.arrow)
                arrows(x0=Idval[iDecision+1],y0=ylim[2],x1=Idval[iDecision],y1=yc[iDecision]+sep.arrow[2],lty=1,length=size.arrow,lwd=lwd.arrow)
            }else{
                arrows(x0=Idval[iDecision],y0=yc[iDecision]-sep.arrow[1],x1=Ival[iDecision+1],y1=ylim[1],lty=1,length=size.arrow,lwd=lwd.arrow)
                arrows(x0=Ival[iDecision+1],y0=ylim[2],x1=Idval[iDecision],y1=yc[iDecision]+sep.arrow[2],lty=1,length=size.arrow,lwd=lwd.arrow)
            }
        }
        arrows(x0=Ival[iDecision+1],y0=ylim[1],x1=Ival[iDecision+1],y1=ylim[2],lty=1,length=size.arrow,lwd=lwd.arrow)
    }
                                        # }}}    
}

