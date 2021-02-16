PlotBoundaries <- function(CalcBndObj,  #object from CalcBoundaries
                           type="Z",    #type of boundaries to plot (Z-statistic (Z), p-value (P), effect (E))
                           Itype="rate"){   #information scale on x-axis (rate or absolute (abs))  
    # {{{ Preliminaries 
    if(Itype=="rate"){
        Ival <- CalcBndObj$Ik/CalcBndObj$Imax
        Idval <- CalcBndObj$Id/CalcBndObj$Imax
        xlb <- "Information fraction"
    } else if(Itype=="abs"){
        Ival <- CalcBndObj$Ik
        Idval <- CalcBndObj$Id
        xlb <- "Information"
    } else {
        stop("Itype not correctly specified")
    }
    if( !(type %in% c("Z","E","P")) )  stop("Incorrect type specified")
    if(type=="P"){
        mydigits <- 3
    }else{
        mydigits=2
    }
    # }}}
    # {{{ Values to plot
    # default (if type="Z")
    xu <- Ival
    xd <- Idval
    yu <- CalcBndObj$uk
    yl <- CalcBndObj$lk
    yc <- CalcBndObj$ck
    yh <- qnorm(1-CalcBndObj$alpha)
    whereleg <- "bottomleft"
    ylab <- "Stopping boundary (Z-statistic)"
    ylim <- c(-1,3)

    if(type=="P"){
        whereleg <- "topleft"
        ylab <- "Stopping boundary (P-value)"
        yu <- 1-pnorm(CalcBndObj$uk)
        yl <- 1-pnorm(CalcBndObj$lk)
        yc <- 1-pnorm(CalcBndObj$ck)
        yh <- CalcBndObj$alpha
        ylim <- c(0,1)    
    }
    if(type=="E"){
        ylab <- "Stopping boundary (Effect estimate)"
        yu <- CalcBndObj$uk/sqrt(CalcBndObj$Ik)
        yl <- CalcBndObj$lk/sqrt(CalcBndObj$Ik)
        yc <- CalcBndObj$c/sqrt(CalcBndObj$Id)
        yh <- qnorm(1-CalcBndObj$alpha)/sqrt(CalcBndObj$Ik[length(CalcBndObj$Ik)])
        ylim <- c(min(yl)-0.5,max(yu)+0.5)
    }
    # }}}
    # {{{ Plot
    plot(xu,yu,type="l",lty=2,lwd=2,ylim=ylim,col="green3",xlab=xlb,ylab=ylab,axes=FALSE,
         xlim=c(0,max(xu)))
    lines(xu,yl,col="red",lwd=2,lty=2)
    points(xd,yc,col="black",pch=19,cex=1.5)
    points(xu,yl,col="red",pch=21,bg="red",cex=1.2)
    points(xu,yu,col="green3",pch=21,bg="green3",cex=1.2)
    #---
    specdec <- function(x,k){ format(round(x,k),nsmall=k)}
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
    legend(whereleg,
           c("Stopping bound for efficacy",
             "Stopping bound for futility",
             "Critical value fixed design",
             "Decision boundary"),
           lty=c(2,2,2,NA),
           pch=c(21,21,NA,19),
           col=c("green3","red","grey","black"),
           pt.bg = c("green3","red",NA,NA),
           bg="white")
    # }}}    
}

#consider whether we want to somehow integrate the plotting function of the Rpact package and add to that, but it
#seems to be a lot of work, since their options for plotting on different scales do not seem to work well
