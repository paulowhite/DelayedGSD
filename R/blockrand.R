### blockrand.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Aug 28 2020 (16:24) 
## Version: 
## Last-Updated: Aug 28 2020 (16:28) 
##           By: Paul Blanche
##     Update #: 7
#----------------------------------------------------------------------
## 
### Commentary: 
##  Block randomization.
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

##' @title Block randomization
##' @param seed
##' @param rand.block
##' @param n
##' @details Nothing
##' @return Vector or binary otcome.
##' @author Paul Blanche
##' 
##' @examples
##' blockrand(123,c(1,1,0,0),10)
blockrand  <-  function(seed,rand.block,n){
    set.seed(seed)
    blocksize <- length(rand.block)
    block  <-  rep(1:ceiling(n/blocksize), each = blocksize)
    a1  <-  data.frame(block, rand=runif(length(block)), envelope= 1: length(block))
    a2  <-  a1[order(a1$block,a1$rand),]
    a2$Z <- rep(rand.block,times = length(block)/blocksize)
    Z  <-  a2[order(a2$envelope),"Z"]
    return(Z[1:n])
}
#----------------------------------------------------------------------
### blockrand.R ends here
