### tests.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Jan  7 2022 (14:36) 
## Version: 
## Last-Updated: Jan  7 2022 (14:46) 
##           By: Paul Blanche
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

#Source Paul's functions
sourceDir("~/research/SeqDesignDelayed/DelayedGSD/R")

bCJ2 <- Method3(rho_alpha=1.345,
                rho_beta=1.345,
                alpha=0.025,
                beta=0.1,
                Kmax=3,
                binding=FALSE,
                Info.max=12,
                InfoR.i=c(3.5,6.75,12)/12,
                InfoR.d=c(5.5,8.75)/12,
                delta=1,  
                abseps = 1e-06, 
                direction="smaller",
                sided=1
                )

bCJ <- NonBindingHJ(rho_alpha=1.345,
                    rho_beta=1.345,
                    alpha=0.025,
                    beta=0.1,
                    Kmax=3,
                    Info.max=12,
                    InfoR.i=c(3.5,6.75,12)/12,
                    InfoR.d=c(5.5,8.75)/12,
                    delta=1,  
                    abseps = 1e-06, 
                    direction="smaller",
                    sided=1
                    )

for (i in 1:length(bCJ)){
    print(all(bCJ[[i]]==bCJ2[[2]]))
}

bCJ2binding <- Method3(rho_alpha=1.345,
                rho_beta=1.345,
                alpha=0.025,
                beta=0.1,
                Kmax=3,
                binding=TRUE,
                Info.max=12,
                InfoR.i=c(3.5,6.75,12)/12,
                InfoR.d=c(5.5,8.75)/12,
                delta=1,  
                abseps = 1e-06, 
                direction="smaller",
                sided=1
                )


#----------------------------------------------------------------------
### tests.R ends here
