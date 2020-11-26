PlotBoundaries <- function(CalcBndObj,  #object from CalcBoundaries
                           type="Z",    #type of boundaries to plot (Z-statistic (Z), p-value (P), effect (E))
                           Itype="rate",
                           main=NULL){   #information scale on x-axis (rate or absolute (abs))
  
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
  
  if(type=="Z"){
    plot(Ival,CalcBndObj$uk,type="l",lty=2,lwd=2,ylim=c(-1,3),col="green",xlab=xlb,ylab="Stopping boundary (Z-statistic)",main=main)
    lines(Ival,CalcBndObj$lk,col="red",lwd=2,lty=2)
    points(Idval,CalcBndObj$ck,col="black",pch=4,cex=1.5)
    points(Ival,CalcBndObj$lk,col="red",pch=21,bg="red",cex=1.2)
    points(Ival,CalcBndObj$uk,col="green",pch=21,bg="green",cex=1.2)
    abline(h=qnorm(1-CalcBndObj$alpha),lty=2,col="grey")
    legend("bottomright",c("Stopping bound for efficacy","Stopping bound for futility","Critical value fixed design","Decision boundary"),lty=c(2,2,2,NA),pch=c(21,21,NA,4),col=c("green","red","grey","black"),pt.bg = c("green","red",NA,NA))
  } else if(type=="P"){
    plot(Ival,1-pnorm(CalcBndObj$uk),type="l",lty=2,lwd=2,ylim=c(0,1),col="green",xlab=xlb,ylab="Stopping boundary (P-value)",main=main)
    lines(Ival,1-pnorm(CalcBndObj$lk),col="red",lwd=2,lty=2)
    points(Idval,1-pnorm(CalcBndObj$ck),col="black",pch=4,cex=1.5)
    points(Ival,1-pnorm(CalcBndObj$lk),col="red",pch=21,bg="red",cex=1.2)
    points(Ival,1-pnorm(CalcBndObj$uk),col="green",pch=21,bg="green",cex=1.2)
    abline(h=CalcBndObj$alpha,lty=2,col="grey")
    legend("topright",c("Stopping bound for efficacy","Stopping bound for futility","Critical value fixed design","Decision boundary"),lty=c(2,2,2,NA),pch=c(21,21,NA,4),col=c("green","red","grey","black"),pt.bg = c("green","red",NA,NA))
  } else if(type=="E"){
    plot(Ival,CalcBndObj$uk/CalcBndObj$Ik,type="l",lty=2,lwd=2,ylim=c(min(CalcBndObj$lk/CalcBndObj$Ik)-0.5,max(CalcBndObj$uk/CalcBndObj$Ik)+0.5),col="green",xlab=xlb,ylab="Stopping boundary (Effect estimate)",main=main)
    lines(Ival,CalcBndObj$lk/CalcBndObj$Ik,col="red",lwd=2,lty=2)
    points(Idval,CalcBndObj$ck/CalcBndObj$Id,col="black",pch=4,cex=1.5)
    points(Ival,CalcBndObj$lk/CalcBndObj$Ik,col="red",pch=21,bg="red",cex=1.2)
    points(Ival,CalcBndObj$uk/CalcBndObj$Ik,col="green",pch=21,bg="green",cex=1.2)
    abline(h=qnorm(1-CalcBndObj$alpha)/CalcBndObj$Ik[length(CalcBndObj$Ik)],lty=2,col="grey")
    legend("bottomright",c("Stopping bound for efficacy","Stopping bound for futility","Critical value fixed design","Decision boundary"),lty=c(2,2,2,NA),pch=c(21,21,NA,4),col=c("green","red","grey","black"),pt.bg = c("green","red",NA,NA))
  } else {
    stop("Incorrect type specified")
  }
  
}

#consider whether we want to somehow integrate the plotting function of the Rpact package and add to that, but it
#seems to be a lot of work, since their options for plotting on different scales do not seem to work well