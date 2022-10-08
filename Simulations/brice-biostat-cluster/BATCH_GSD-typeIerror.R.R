## * Header 
## path <- "x:/DelayedGSD/"
## setwd(path)
## source("BATCH_GSD-typeIerror.R.R")
## sbatch -a 1-1 -J 'GSD-typeIerror.R' --output=/dev/null --error=/dev/null R CMD BATCH --vanilla BATCH_GSD-typeIerror.R.R /dev/null 

rm(list = ls())
gc()

## * seed
iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
if(is.na(iter_sim)){iter_sim <- 1}
if(is.na(n.iter_sim)){n.iter_sim <- 10}

## * prepare export
path <- "."
path.res <- file.path(path,"Results","GSD-typeIerror.R")
if(dir.exists(path.res)==FALSE){
    if(dir.exists(file.path(path,"Results"))==FALSE){
    dir.create(file.path(path,"Results"))
    }
    dir.create(path.res)
}
path.output <- file.path(path,"output","GSD-typeIerror.R")
if(dir.exists(path.output)==FALSE){
    if(dir.exists(file.path(path,"output"))==FALSE){
    dir.create(file.path(path,"output"))
    }
    dir.create(path.output)
}

## * libraries
library(DelayedGSD)
source("FCT.R")

## * settings
n.sim <- 50

## * function to execute
## method2num[c(4,10,14),]
##    index method binding correction fixC
## 4      4      1   FALSE      FALSE TRUE
## 10    10      2   FALSE      FALSE TRUE
## 14    14      3   FALSE      FALSE TRUE


out <- FCT_simGSD(method = 4,## c(4,10,14),
                  kMax = 2,
                  InfoR.i = c(0.5,1),
                  InfoR.d = c(0.55,1),
                  PropForInterim = 0.5,
                  deltaPlan = c(0,0.6,0.8),
                  deltaTrue = c(0,0,0),
                  allsd = c(2.5,2.1,2.4),
                  mean0 = c(10,0,0),
                  n.sample = 1000,
                  path = path.res,
                  export.tempo = TRUE,
                  n.sim = c(n.sim,iter_sim),
                  seed = 140786598)


## * R version
print(sessionInfo())

## * gather and process results
if(FALSE){
path <- "x:/DelayedGSD/"
setwd(path)

path.GSD-typeIerror.R <- file.path("Results","GSD-typeIerror.R")
allRes.tempo <- butils::sinkDirectory(path.GSD-typeIerror.R, string.keep = "tempo")
allRes.final <- butils::sinkDirectory(path.GSD-typeIerror.R, string.exclude = "tempo")
}

	
