### PointwiseEstimate.R --- 
#----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Jul  6 2021 (10:16) 
## Version: 
## Last-Updated: Jul  6 2021 (16:21) 
##           By: Paul Blanche
##     Update #: 244
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

rm(list=ls())
# Note: delta is the treatment effect on the Z scale, i.e. assuming sd=1


pk <- function(sk,k,AllIk,delta,Alluk,Alllk,K,TRACE=FALSE){
    if(TRACE) print("start pk")
    # {{{ check
    ## length of Allsk, AllIk, Alluk, Alllk should be k !!!
    # }}}
    
    Ik <- AllIk[k]

    ## sk <- Allsk[k]
    if(k>1){
        ## sM1 <- Allsk[k-1]
        IkM1 <- AllIk[k-1]
        lkM1 <- Alllk[k-1]
        ukM1 <- Alluk[k-1]
    }

    
    # {{{ k=1
    if(k==1){
        out <- dnorm(sk-delta*sqrt(Ik)) # to double check
    }
    # }}}

    # {{{ k
    if(k>1){
        ## tauk <- Ik-IkM1
        ToInt <- function(s){
            ## NewAllsk <-  c(Allsk[1:(k-2)],s)
            ## ss <- ((sk-s)-delta*(sqrt(Ik)-sqrt(IkM1)))/sqrt(Ik)
            ss <- ((sk-s)-delta*(sqrt(Ik)-sqrt(IkM1)))/sqrt(Ik-IkM1)
            
            ss <- ((sk-s)-delta*(Ik-IkM1)/I1)/sqrt(Ik-IkM1)
            ## pk(s,k-1,AllIk[1:(k-1)],delta,Alluk[1:(k-1)],Alllk[1:(k-1)],K)*(1/sqrt(Ik))*dnorm(ss)
            pk(s,k-1,AllIk[1:(k-1)],delta,Alluk[1:(k-1)],Alllk[1:(k-1)],K)*dnorm(ss)*sqrt(Ik/(Ik-IkM1))
        }
        ## out <- integrate(ToInt,lower=-Inf,upper=lkM1)$value + integrate(ToInt,lower=ukM1,upper=Inf)$value           
        out <- integrate(ToInt,lower=lkM1,upper=ukM1)$value           
    }
    # }}}    
    out
}

fk <- function(sk,k,AllIk,delta,Alluk,Alllk,K,TRACE=FALSE){
    if(TRACE) print("start fk")
    # {{{ check
    ## length of Allsk, AllIk, Alluk, Alllk should be k !!!
    # }}}

    # {{{ parameters
    ## Ik <- AllIk[k]
    ## sk <- Allsk[k]
    lk <- Alllk[k]
    uk <- Alluk[k]
    # }}}   
    out <- 0
    ## if(k==2 )
    ## browser()
    if(sk >uk | sk < lk | k==K){
        out <- pk(sk,k,AllIk,delta,Alluk,Alllk,K)
    }
    out
}


#------- Check densities ----------
MyDelta <- 0.9915
MyAllIk <- c(1:4)*22
MyAlllk <- -4.049/sqrt(1:4)
MyAlluk <- 4.049/sqrt(1:4)

xxx <- seq(from=-20,to=20,length.out=5000)
## yyy1 <- sapply(xxx,function(s) fk(s,k=1,MyAllIk[1],delta=MyDelta,MyAlluk[1],MyAlllk[1],K=4,TRACE=FALSE))
## yyy2 <- sapply(xxx,function(s) fk(s,k=2,MyAllIk[1:2],delta=MyDelta,MyAlluk[1:2],MyAlllk[1:2],K=4,TRACE=FALSE))
## yyy3 <- sapply(xxx,function(s) fk(s,k=3,MyAllIk[1:3],delta=MyDelta,MyAlluk[1:3],MyAlllk[1:3],K=4,TRACE=FALSE))
## yyy4 <- sapply(xxx,function(s) fk(s,k=4,MyAllIk[1:4],delta=MyDelta,MyAlluk[1:4],MyAlllk[1:4],K=4,TRACE=FALSE))

yyy1 <- sapply(xxx,function(s) fk(s,k=1,MyAllIk[1],delta=MyDelta,MyAlluk[1],MyAlllk[1],K=4,TRACE=FALSE))
yyy2 <- sapply(xxx,function(s) fk(s,k=2,MyAllIk[1:2],delta=MyDelta,MyAlluk[1:2],MyAlllk[1:2],K=4,TRACE=FALSE))
yyy3 <- sapply(xxx,function(s) fk(s,k=3,MyAllIk[1:3],delta=MyDelta,MyAlluk[1:3],MyAlllk[1:3],K=4,TRACE=FALSE))
yyy4 <- sapply(xxx,function(s) fk(s,k=4,MyAllIk[1:4],delta=MyDelta,MyAlluk[1:4],MyAlllk[1:4],K=4,TRACE=FALSE))


plot(xxx,yyy1,col=1,type="l",lwd=2,ylim=range(c(yyy1,yyy2,yyy3,yyy4)))
abline(v=MyDelta,lty=3)
## abline(v=MyDelta*MyAllIk[1],lty=3)
lines(xxx,yyy2,col=2,lwd=2)
## abline(v=MyDelta*MyAllIk[2],lty=3)
lines(xxx,yyy3,col=3,lwd=2)
## abline(v=MyDelta*MyAllIk[3],lty=3)
lines(xxx,yyy4,col=4,lwd=2)
## abline(v=MyDelta*MyAllIk[4],lty=3)
legend("top",col=1:4,legend=1:4,lwd=2)

f1 <- function(s) fk(s,k=1,MyAllIk[1:1],delta=MyDelta,MyAlluk[1:1],MyAlllk[1:1],K=4,TRACE=FALSE)
f2 <- function(s) fk(s,k=2,MyAllIk[1:2],delta=MyDelta,MyAlluk[1:2],MyAlllk[1:2],K=4,TRACE=FALSE)
f3 <- function(s) fk(s,k=3,MyAllIk[1:3],delta=MyDelta,MyAlluk[1:3],MyAlllk[1:3],K=4,TRACE=FALSE)
f4 <- function(s) fk(s,k=4,MyAllIk[1:4],delta=MyDelta,MyAlluk[1:4],MyAlllk[1:4],K=4,TRACE=FALSE)

v.f1 <- Vectorize(f1)
v.f2 <- Vectorize(f2)
v.f3 <- Vectorize(f3)
v.f4 <- Vectorize(f4)

int1 <- integrate(v.f1,lower=-Inf,upper=Inf)$value
int2 <- integrate(v.f2,lower=-Inf,upper=Inf)$value
int3 <- integrate(v.f3,lower=-Inf,upper=Inf)$value
int4 <- integrate(v.f4,lower=-Inf,upper=Inf)$value

Allint <- c(int1,int2,int3,int4)
sum(Allint)


#--------------

BiasML <- function(AllIk,delta,Alluk,Alllk,K,TRACE=FALSE){
    ToInt <- function(s,k) (s/AllIk[k])*fk(s,k,AllIk[1:k],delta,Alluk[1:k],Alllk[1:k],K,TRACE)
    Bias <- 0
    # Step k=1,...,K-1
    browser()
    for(k in 1:(K-1)){
        if(TRACE) print(paste0("k=",k," in Bias computation"))
        ToIntk <- function(s) ToInt(s,k)
        v.ToIntk <- Vectorize(ToIntk)
        ## browser()
        integrate(v.ToIntk,lower=-5,upper=5)$value
        
        Bias <- Bias +  integrate(v.ToIntk,lower=-Inf,upper=Alllk[k])$value +  integrate(v.ToIntk,lower=Alluk[k],upper=Inf)$value    
    }
    # Step k=K
    ## browser()
    if(TRACE) print(paste0("k=",K," in Bias computation"))
    ToIntK <- function(s) (s/AllIk[K])*fk(s,k,AllIk[1:K],delta,Alluk[1:K],Alllk[1:K],K,TRACE)
    v.ToIntK <- Vectorize(ToIntK)
    Bias <- Bias + integrate(v.ToIntK,lower=-Inf,upper=Inf)$value
    Bias - delta
}


# Example book page 96, fig 4.2

BiasML(AllIk=MyAlllk,
       delta=MyDelta,
       Alluk=MyAlluk,
       Alllk=MyAlllk,
       K=4,TRACE=TRUE)



# To do
# Check density condition as p 23



#----------------------------------------------------------------------
### PointwiseEstimate.R ends here
