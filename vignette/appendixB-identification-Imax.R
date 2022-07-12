### appendixB-identification-Imax.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul 12 2022 (09:12) 
## Version: 
## Last-Updated: jul 12 2022 (18:48) 
##           By: Brice Ozenne
##     Update #: 20
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(truncnorm)
library(numDeriv)
library(DelayedGSD)
library(tmvnorm)

## * Setting
## CalcBoundaries(kMax = 2)

l1 <- 0.3788015
u1 <- 2.4977055
u2 <- 2.2
l2 <- 1.1

n1 <- 50
n2 <- 50
w1 <- sqrt(n1/(n1+n2))
w2 <- sqrt(n2/(n1+n2))

Sigma <- matrix(c(1,w1,w1,1),2,2)
alpha2 <- 0.02500 - 0.00625
beta2 <- 0.20 - 0.05

Imax <- 3.64016
delta <- 0.15

## * Lemma B.1

## ** Value
GS.B1 <- pmvnorm(lower = c(l1,u2),
                 upper = c(u1,Inf),
                 mean = rep(0,2),
                 sigma= Sigma)
GS.B1
## [1] 0.01111284

integrate(f = function(x){
    pnorm((u2 - w1*x)/w2,lower.tail = FALSE)*dnorm(x)
}, lower = l1, upper = u1, rel.tol = 1e-8)$value
## [1] 0.01111284

(pnorm(u1)-pnorm(l1))-integrate(f = function(x){pnorm((u2 - w1*x)/w2)*dnorm(x)}, lower = l1, upper = u1, rel.tol = 1e-8)$value
## [1] 0.01111284

(pnorm(u1)-pnorm(l1))*(1-integrate(f = function(x){pnorm((u2 - w1*x)/w2)*dtruncnorm(x, a = l1, b = u1)}, lower = -Inf, upper = Inf, rel.tol = 1e-8)$value)
## [1] 0.01111284

FUN.I <- function(x,y){
    (w2*x+w1*y<=u2)*dnorm(x)*dtruncnorm(y, a = l1, b = u1)
}
(pnorm(u1)-pnorm(l1))*(1-integrate(f = function(x){sapply(x, FUN = function(iX){integrate(f = FUN.I, x = iX, lower = l1, upper = u1, rel.tol = 1e-7)$value})}, lower = -10, upper = 10, rel.tol = 1e-7)$value)
## [1] 0.01111283

library(pracma)
(pnorm(u1)-pnorm(l1))*(1-integral2(fun = FUN.I, xmin = -10, xmax = +10, ymin = l1, ymax = u1, reltol = 1e-8, maxlist = 1e4)$Q)
## [1] 0.0101619

library(cubature)
FUN.I2 <- function(x){
    (w2*x[1]+w1*x[2]<=u2)*dnorm(x[1])*dtruncnorm(x[2], a = l1, b = u1)
}
res <- (pnorm(u1)-pnorm(l1))*(1-hcubature(f = FUN.I2, lowerLimit = c(-10,l1), upperLimit = c(10,u1), tol = 1e-6)$integral)
res
## [1] 0.01111361

## ** Derivative
Psi <- function(ZZ){
    integrate(f = function(x){pnorm((ZZ - w1*x)/w2)*dtruncnorm(x, a = l1, b = u1)}, lower = -Inf, upper = Inf, rel.tol = 1e-8)$value
}

findU <- function(L){uniroot(f = function(X){
    pmvnorm(lower = c(L,X),
            upper = c(2.4977055,Inf),
            mean = rep(0,2),
            sigma= Sigma)-alpha2
}, lower = 0, upper = 3)$root}
findU(l1)
## [1] 2.000239

findU2 <- function(L){uniroot(f = function(X){
    (pnorm(u1)-pnorm(l1))*(1-Psi(L))-alpha2
}, lower = 0, upper = 3)$root}
findU2(l1)

jacobian(findU, l1)
##             [,1]
## [1,] -0.05680709

-alpha2*dnorm(l1)/(pnorm(u1)-pnorm(l1))^2



jacobian(Psi, u2)
##            [,1]
## [1,] 0.08819166
dPsi <- function(ZZ){
    integrate(f = function(x){dnorm((ZZ - w1*x)/w2)*dtruncnorm(x, a = l1, b = u1)}, lower = -Inf, upper = Inf, rel.tol = 1e-8)$value
}
1/w2*dPsi(u2)
## [1] 0.08819166
dPsi <- function(ZZ){
    integrate(f = function(x){dnorm((ZZ - w1*x)/w2)*dnorm(x)}, lower = l1, upper = u1, rel.tol = 1e-8)$value
}
1/(w2*(pnorm(u1)-pnorm(l1)))*dPsi(u2)
## [1] 0.08819166
dPsi <- function(ZZ){
    integrate(f = function(x){dnorm((ZZ - w1*x)/w2)*dnorm(x)}, lower = l1, upper = u1, rel.tol = 1e-8)$value
}
1/(w2*(pnorm(u1)-pnorm(l1)))*dPsi(u2)
## [1] 0.08819166
dPsi <- function(ZZ){
    integrate(f = function(x){dnorm(x*sqrt(1+w1^2/w2^2)-w1/w2*ZZ/(w2*sqrt(1+w1^2/w2^2)))*exp(-0.5*ZZ^2/(w2^2*(1+w1^2/w2^2)))/sqrt(2*pi)}, lower = l1, upper = u1, rel.tol = 1e-8)$value
}
1/(w2*(pnorm(u1)-pnorm(l1)))*dPsi(u2)
## [1] 0.08819166
mu <- w1*u2/(w2^2+w1^2)
sigma <- 1/sqrt(1+w1^2/w2^2)
1/(w2*sqrt(2*pi)*(pnorm(u1)-pnorm(l1)))*(pnorm(u1, mean = mu, sd = sigma)-pnorm(l1, mean = mu, sd = sigma))*exp(-0.5*u2^2/(w2^2*(1+w1^2/w2^2)))*sigma
1/(w2*sqrt(2*pi)*(pnorm(u1)-pnorm(l1)))*(pnorm(u1, mean = mu, sd = sigma)-pnorm(l1, mean = mu, sd = sigma))*exp(-0.5*u2^2/(w2^2*(1+w1^2/w2^2)))*sigma
1/(w2*(pnorm(u1)-pnorm(l1)))*(pnorm(u1, mean = mu, sd = sigma)-pnorm(l1, mean = mu, sd = sigma))*dnorm(u2, mean = 0, sd = sigma)


## * Lemma B.2

## ** Value
GS.B2 <- pmvnorm(lower = c(l1,-Inf),
                 upper = c(u1,l2),
                 mean = c(delta*sqrt(n1),delta*(w1*sqrt(n1)+w2*sqrt(n2))),
                 sigma = Sigma)
GS.B2
## [1] 0.1614915

FUN.Ibis <- function(x){
    (l1<=x[1])*(x[1]<=u1)*(x[2]<=l2)*dmvnorm(cbind(x[1],x[2]),
                                             mean = c(delta*sqrt(n1),w1*delta*sqrt(n1)+w2*delta*sqrt(n2)),
                                             sigma = Sigma)
}
res <- cubature::hcubature(f = FUN.Ibis, lowerLimit = c(l1,-10), upperLimit = c(u1,10))$integral
res
## > [1] 0.1614915

FUN.Ibis <- function(x){
    (l1<=x[1])*(x[1]<=u1)*(w1*x[1]+w2*x[2]<=l2)*dnorm(x[1], mean = delta*sqrt(n1))*dnorm(x[2], mean = delta*sqrt(n2))
}
FUN.Ibis <- function(x){
    (x[2]<=(l2-w1*x[1])/w2)*dnorm(x[1], mean = delta*sqrt(n1))*dnorm(x[2], mean = delta*sqrt(n2))
}
res <- cubature::hcubature(f = FUN.Ibis, lowerLimit = c(l1,-10), upperLimit = c(u1,10))$integral
res
## [1] 0.1608881

integrate(f = function(x){
    pnorm((l2-w1*x)/w2, mean = delta*sqrt(n2))*dnorm(x, mean = delta*sqrt(n1))
}, lower = l1, upper = u1, rel.tol = 1e-8)$value
## [1] 0.1614915

integrate(f = function(x){
    pnorm((l2 - w2*delta*sqrt(n2) - w1*x) / w2)*dnorm(x - delta*sqrt(n1))
}, lower = l1, upper = u1, rel.tol = 1e-8)$value
## [1] 0.1614915
integrate(f = function(x){
    pnorm((l2 - w2*delta*sqrt(n2) - w1*delta*sqrt(n1) - w1*x) / w2)*dnorm(x)
}, lower = l1-delta*sqrt(n1), upper = u1-delta*sqrt(n1), rel.tol = 1e-8)$value
## [1] 0.1614915
(pnorm(u1-delta*sqrt(n1))-pnorm(l1-delta*sqrt(n1)))*integrate(f = function(x){
    pnorm((l2 - w2*delta*sqrt(n2) - w1*delta*sqrt(n1) - w1*x) / w2)*dtruncnorm(x, a = l1-delta*sqrt(n1), b = u1-delta*sqrt(n1))
}, lower = -Inf, upper = Inf, rel.tol = 1e-8)$value
## [1] 0.1614915

FUN.I <- function(x,y){
    (w2*(x+delta*sqrt(n2))+w1*(y+delta*sqrt(n1))<=l2)*dnorm(x)*dtruncnorm(y, a = l1-delta*sqrt(n1), b = u1-delta*sqrt(n1))
}
res <- (pnorm(u1-delta*sqrt(n1))-pnorm(l1-delta*sqrt(n1)))*integrate(f = function(x){sapply(x, FUN = function(iX){integrate(f = FUN.I, x = iX, lower = l1-delta*sqrt(n1), upper = u1-delta*sqrt(n1), rel.tol = 1e-7)$value})}, lower = -10, upper = 10, rel.tol = 1e-7)$value
res
## [1] 0.1608881

##----------------------------------------------------------------------
### appendixB-identification-Imax.R ends here
