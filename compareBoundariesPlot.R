#Function to automatically load all R functions in a directory
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

#Source all functions
sourceDir("R")

#Method 1
M1nofixCbinding <-  CalcBoundaries(kMax=3,
                                 InfoR.i=c(0.5,0.75,1),
                                 method=1,
                                 cNotBelowFixedc = F,
                                 InfoR.d = c(0.6,0.85),
                                 bindingFutility = T)

M1fixCbinding <-  CalcBoundaries(kMax=3,
                                   InfoR.i=c(0.5,0.75,1),
                                   method=1,
                                   cNotBelowFixedc = T,
                                   InfoR.d = c(0.6,0.85),
                                   bindingFutility = T)

M1fixCnonbinding <-  CalcBoundaries(kMax=3,
                                 InfoR.i=c(0.5,0.75,1),
                                 method=1,
                                 cNotBelowFixedc = T,
                                 InfoR.d = c(0.6,0.85),
                                 bindingFutility = F)

M1nofixCnonbinding <-  CalcBoundaries(kMax=3,
                                 InfoR.i=c(0.5,0.75,1),
                                 method=1,
                                 cNotBelowFixedc = F,
                                 InfoR.d = c(0.6,0.85),
                                 bindingFutility = F)

#Method 2
M2nofixCbinding <-  CalcBoundaries(kMax=3,
                                   InfoR.i=c(0.5,0.75,1),
                                   method=2,
                                   cNotBelowFixedc = F,
                                   InfoR.d = c(0.6,0.85),
                                   bindingFutility = T)

M2fixCbinding <-  CalcBoundaries(kMax=3,
                                 InfoR.i=c(0.5,0.75,1),
                                 method=2,
                                 cNotBelowFixedc = T,
                                 InfoR.d = c(0.6,0.85),
                                 bindingFutility = T)

M2fixCnonbinding <-  CalcBoundaries(kMax=3,
                                    InfoR.i=c(0.5,0.75,1),
                                    method=2,
                                    cNotBelowFixedc = T,
                                    InfoR.d = c(0.6,0.85),
                                    bindingFutility = F)

M2nofixCnonbinding <-  CalcBoundaries(kMax=3,
                                      InfoR.i=c(0.5,0.75,1),
                                      method=2,
                                      cNotBelowFixedc = F,
                                      InfoR.d = c(0.6,0.85),
                                      bindingFutility = F)

#Method 3
M3binding <-  CalcBoundaries(kMax=3,
                                 InfoR.i=c(0.5,0.75,1),
                                 method=3,
                                 cNotBelowFixedc = T,
                                 InfoR.d = c(0.6,0.85),
                                 bindingFutility = T)

M3nonbinding <-  CalcBoundaries(kMax=3,
                                    InfoR.i=c(0.5,0.75,1),
                                    method=3,
                                    cNotBelowFixedc = T,
                                    InfoR.d = c(0.6,0.85),
                                    bindingFutility = F)

#comparing all options Method 1
par(mfrow=c(2,2))
plot(M1fixCbinding,legend=F)
plot(M1nofixCbinding,legend=F)
plot(M1fixCnonbinding,legend=F)
plot(M1nofixCnonbinding,legend=F)

#comparing all options Method 2
par(mfrow=c(2,2))
plot(M2fixCbinding,legend=F)
plot(M2nofixCbinding,legend=F)
plot(M2fixCnonbinding,legend=F)
plot(M2nofixCnonbinding,legend=F)

#comparing all options Method 3
par(mfrow=c(1,2))
plot(M3binding,legend=F)
plot(M3nonbinding,legend=F)

#comparing Method 1,2,3 binding, free choice of c (when allowed)
par(mfrow=c(1,3))
plot(M1nofixCbinding,legend=F)
plot(M2nofixCbinding,legend=F)
plot(M3binding,legend=F)

#comparing Method 1,2,3 binding, fixed c
par(mfrow=c(1,3))
plot(M1fixCbinding,legend=F)
plot(M2fixCbinding,legend=F)
plot(M3binding,legend=F)

#comparing Method 1,2,3 non-binding, free choice of c (when allowed)
par(mfrow=c(1,3))
plot(M1nofixCnonbinding,legend=F)
plot(M2nofixCnonbinding,legend=F)
plot(M3nonbinding,legend=F)

#comparing Method 1,2,3 non-binding, fixed c
par(mfrow=c(1,3))
plot(M1fixCnonbinding,legend=F)
plot(M2fixCnonbinding,legend=F)
plot(M3nonbinding,legend=F)
