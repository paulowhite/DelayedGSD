### 0-onload.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 17 2022 (16:45) 
## Version: 
## Last-Updated: feb 17 2022 (16:46) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

.onAttach <- function(lib, pkg="DelayedGSD") {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}


##----------------------------------------------------------------------
### 0-onload.R ends here
