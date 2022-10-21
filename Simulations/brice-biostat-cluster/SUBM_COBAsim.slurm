#!/bin/bash

#SBATCH --job-name=COBAsim  # Job name
#SBATCH --array=1-100     # Number of jobs
#SBATCH --partition=standard    # Queue. Use long for run time >8h and standard otherwise
#SBATCH --mem=1G # Memory limit, e.g. reserve 1 GB memory 
#SBATCH --output=output/slurm/slurm-%x-%a-%J.out # Direct output to subdirectory
#SBATCH --error=output/slurm/slurm-%x-%a-%J.out # Same file for output and error

R CMD BATCH --vanilla BATCH_missing_power.R output/missing_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla BATCH_missing_typeI.R output/missing_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla BATCH_noFixC_power.R output/noFixC_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla BATCH_noFixC_typeI.R output/noFixC_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla BATCH_nomissing_power.R output/nomissing_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla BATCH_nomissing_typeI.R output/nomissing_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla BATCH_nonbinding_power.R output/nonbinding_power/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout
R CMD BATCH --vanilla BATCH_nonbinding_typeI.R output/nonbinding_typeI/R-$SLURM_JOB_NAME-$SLURM_ARRAY_TASK_ID-$SLURM_JOB_ID.Rout

## go to directory    ## cd c:/Users/hpl802/Documents/Github/DelayedGSD/Simulations/brice-biostat-cluster/
                      ## cd /projects/biostat01/people/hpl802/DelayedGSD/
## clean outputs      ## rm -r output/slurm/*
                      ## rm -r output/missing_power/*
                      ## rm -r output/missing_typeI/*
                      ## rm -r output/noFixC_power/*
                      ## rm -r output/noFixC_typeI/*
                      ## rm -r output/nomissing_power/*
                      ## rm -r output/nomissing_typeI/*
                      ## rm -r output/nonbinding_power/*
                      ## rm -r output/nonbinding_typeI/*
## clean results      ## rm -r Results/missing_power/*
                      ## rm -r Results/missing_typeI/*
                      ## rm -r Results/noFixC_power/*
                      ## rm -r Results/noFixC_typeI/*
                      ## rm -r Results/nomissing_power/*
                      ## rm -r Results/nomissing_typeI/*
                      ## rm -r Results/nonbinding_power/*
                      ## rm -r Results/nonbinding_typeI/*
## submission command ## sbatch SUBM_COBAsim.slurm

## submission output  ## Submitted batch job 2388429 (time using C-c t: torsdag 22-10-20 at 15:59) 

## init               ## module module load gcc/11.2.0 R/4.1.2
                      ## cat .bash_profile (in the root directory on the server to visualize the current configuration)
                      ## vi .bash_profile (in the root directory on the server to modify the current configuration)

## commands           ## squeue           : view current jobs 
                      ## squeue -u id     : view current jobs for user id
                      ## scancel          : delete the job with job id 1034 from the queue type
                      ## sinfo            : view state of the queues
                      ## sshare           : check own usage of the server
                        
## documentation      ## from SLURM: https://slurm.schedmd.com/pdfs/summary.pdf: key functions
                      ##           : https://slurm.schedmd.com/sbatch.html (filename pattern)
                      ## from KU   : https://hpc.ku.dk/documentation/slurm.html
                      ##           : https://github.com/aejensen/computing-biostat-ucph/blob/main/README.md
