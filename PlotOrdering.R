rm(list=ls())

if(system("whoami",intern=TRUE)=="paul"){  
    setwd("~/research/SeqDesignDelayed/DelayedGSD/")
}else if(system("whoami",intern=TRUE)=="brice"){  
    setwd("~/Documents/GitHub/DelayedGSD/")
}
sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
        if(trace) cat(nm,":")
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}
sourceDir("Rfunctions/")


bnds <- CalcBoundaries(kMax=3,informationRates = c(0.5,0.75,1),Id=c(0.55,0.8))

plot(bnds$Ik,bnds$lk,col="red",pch=21,bg="red",ylim=c(0,max(bnds$uk+0.5)),xlab="Information",ylab="Stopping boundary (Z-statistic)")
lines(bnds$Ik,bnds$lk,col="red",lty=2)
points(bnds$Ik,bnds$uk,col="green",lty=2,pch=21,bg="green")
lines(bnds$Ik,bnds$uk,col="green",lty=2)
points(bnds$Id,bnds$ck,pch=4)
arrows(x0=bnds$Id[1],y0=0,x1=bnds$Id[1],y1=bnds$ck[1]-0.1,lty=1,length=0.1)
arrows(x0=bnds$Id[1],y0=bnds$ck[1]-0.1,x1=bnds$Id[2],y1=0,lty=1,length=0.1)
arrows(x0=bnds$Id[2],y0=0,x1=bnds$Id[2],y1=bnds$ck[2]-0.1,lty=1,length=0.1)
arrows(x0=bnds$Id[2],y0=bnds$ck[2]-0.1,x1=bnds$Ik[3],y1=0,lty=1,length=0.1)
arrows(x0=bnds$Ik[3],y0=0,x1=bnds$Ik[3],y1=3,lty=1,length=0.1)
arrows(x0=bnds$Ik[3],y0=3,x1=bnds$Id[2],y1=bnds$ck[2]+0.1,lty=1,length=0.1)
arrows(x0=bnds$Id[2],y0=bnds$ck[2]+0.1,x1=bnds$Id[2],y1=3,lty=1,length=0.1)
arrows(x0=bnds$Id[2],y0=3,x1=bnds$Id[1],y1=bnds$ck[1]+0.1,lty=1,length=0.1)
arrows(x0=bnds$Id[1],y0=bnds$ck[1]+0.1,x1=bnds$Id[1],y1=3,lty=1,length=0.1)
