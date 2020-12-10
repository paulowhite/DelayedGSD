## * R packages
## https://cran.r-project.org/web/views/NumericalMathematics.html
## Section: Root Finding and Fixed Points
library(BB)
## For solving nonlinear systems of equations the BB package provides Barzilai-Borwein spectral methods in sane(), including a derivative-free variant in dfsane(), and multi-start features with sensitivity analysis.
library(nleqslv)
##Package nleqslv solves nonlinear systems of equations using alternatively the Broyden or Newton method, supported by strategies such as line searches or trust regions. 

## * Solvers
##' @title One dimensional root finding
##' @description Find the zero of a function by specifying a intial guess (dfsane, L-BFGS-B, nleqslv) or a range of possible values (uniroot) or both (L-BFGS-B).
##' All algorithms but nleqslv are derivative free.
##'
##' @param FUN [function] univariate function for which the 0 should be located.
##' @param method [character] algorithm to be used: one of \code{"dfsane"}, \code{"L-BFGS-B"}, \code{"nleqslv"}, \code{"uniroot"}.
##' @param start [numeric] starting value.
##' @param lower [numeric] lower possible value.
##' @param upper [numeric] largest possible value.
##' @param maxiter [integer] maximum number of iterations.
##' @param tol [numeric] acceptable error in FUN when locating the 0.
rootSolve <- function(FUN, method, n.try = 1,
                      start = NULL, lower = NULL, upper = NULL,
                      maxiter = 1000, tol = 1e-7){

    if(n.try > 1){       
        res <- do.call(rbind,lapply(1:n.try, function(iTry){
            rootSolve(FUN = FUN, method = method, n.try = 1, start = start, lower = lower, upper = upper)
        }))
        out <- res[which.min(res$abs.error),]
        attr(out,"all") <- res
    }else{

        method <- match.arg(method, c("uniroot","L-BFGS-B","dfsane","nleqslv"))
        if(method == "uniroot"){
            if(is.null(lower) || is.null(upper)){
                stop("Arguments \'lower\' and \'upper\' must be specified when using the uniroot method. \n")
            }
            if(!is.null(start)){
                stop("Arguments \'start\' must NOT be specified when using the uniroot method. \n")
            }
            res <- stats::uniroot(f = FUN, lower = lower, upper = upper, maxiter = maxiter, tol = tol)
            out <- data.frame(solution = res$root,
                              abs.error = abs(res$f.root),
                              convergence = !is.na(res$estim.prec),
                              method = method
                              )
        }else if(method == "dfsane"){
            if(is.null(start)){
                stop("Arguments \'start\' must be specified when using the dfsane method. \n")
            }
            if(!is.null(lower) || !is.null(upper)){
                stop("Arguments \'lower\' and \'upper\' must NOT be specified when using the dfsane method. \n")
            }
            res <- BB::dfsane(fn = FUN, par = start, control = list(maxit = maxiter, tol = tol), quiet = TRUE)
            out <- data.frame(solution = res$par,
                              abs.error = abs(res$residual),
                              convergence = res$convergence==0,
                              method = method
                              )
        }else if(method == "nleqslv"){
            if(is.null(start)){
                stop("Arguments \'start\' must be specified when using the nleqslv method. \n")
            }
            if(!is.null(lower) || !is.null(upper)){
                stop("Arguments \'lower\' and \'upper\' must NOT be specified when using the nleqslv method. \n")
            }
            res <- nleqslv::nleqslv(fn = FUN, x = start, control = list(maxit = maxiter, ftol = tol))
            out <- data.frame(solution = res$x,
                              abs.error = abs(res$fvec),
                              convergence = res$termcd==1,
                              method = method
                              )
        }else if(method == "L-BFGS-B"){
            if(is.null(start)){
                stop("Arguments \'start\' must be specified when using the L-BFGS-B method. \n")
            }
            if(is.null(lower)){lower <- -Inf}
            if(is.null(upper)){upper <- Inf}
            res <- stats::optim(fn = function(x){abs(FUN(x))}, par = start, lower = lower, upper = upper, 
                                method = "L-BFGS-B", control = list(maxit = maxiter, factr = tol))
            out <- data.frame(solution = res$par,
                              abs.error = abs(res$value),
                              convergence = res$convergence==0,
                              method = method
                              )
        }
    }

        return(out)
}


## * Example 1
p <- 5
Sigma <- matrix(apply(expand.grid(1:(p+1),1:(p+1)), 1, function(x){min(x)/max(x)}), nrow = p+1, ncol = p+1)
l <- rep(-2,p)
u <- rep(2,p)


calcP <- function(x){
    if(length(x)>1){
        return(sapply(x,calcP))
    }else{
        a <- pmvnorm(lower = c(l[0:(p-1)],u[p],-Inf),
                     upper = c(u[0:(p-1)],Inf,x),
                     mean=rep(0,p+1),
                     sigma= Sigma) 
        b <- pmvnorm(lower = c(l[0:(p-1)],-Inf,x),
                     upper = c(u[0:(p-1)],l[p],Inf),
                     mean=rep(0,p+1),
                     sigma= Sigma)
        return(as.double(a-b))
    }    
}

## ** solve
## Using different methods
rootSolve(calcP, start = 1, method = "dfsane", n.try = 5)
rootSolve(calcP, lower = -2, upper = 2, method = "uniroot")
rootSolve(calcP, start = 1, method = "L-BFGS-B")
rootSolve(calcP, start = 1, lower = -2, upper = 2, method = "L-BFGS-B")
rootSolve(calcP, start = 1, method = "nleqslv", n.try = 5)

n.try <- 10
ls.tps <- list()
ls.res <- list()
ls.tps$dfsane <- system.time({ls.res$dfsane <- attr(rootSolve(calcP, start = 1, method = "dfsane", n.try = n.try),"all")})
ls.tps$uniroot <- system.time({ls.res$uniroot <- attr(rootSolve(calcP, lower = -2, upper = 2, method = "uniroot", n.try = n.try),"all")})
ls.tps$BFGS <- system.time({ls.res$BFGS <- attr(rootSolve(calcP, start = 1, method = "L-BFGS-B",n.try),"all")})
## ls.tps$L_BFGS_B <- system.time({L_BFGS_B <- rootSolve(calcP, start = 1, lower = -2, upper = -2, method = "L-BFGS-B",n.try)})
ls.tps$nlqslv <- system.time({ls.res$nlqslv <- attr(rootSolve(calcP, start = 1, method = "nleqslv", n.try = n.try),"all")})

M.tps <- do.call(rbind,ls.tps)
M.res <- do.call(rbind,ls.res)
boxplot(solution ~ method, data = M.res)
points(solution ~ method, data = M.res)
