#
# This calls all scripts in the right order to replicate results from my paper
# 'Competing Constructions for German Measure Phrases' (working title)
#
# Warning: This takes quite a long times to run, especially with 
#  =>  model coef & LR BOOTSTRAPPING
#  =>  BAYESIAN / MCMC estimation

# Comment this line if you want to keep your old workspace.
rm(list = ls())

# Set this to the location of R script sources.
setwd('/Users/user/Workingcopies/measurenps/GitHub/R')

# OPTIONS

set.seed(2398651)            # Do NOT change, or results will differ.
save.persistent <- F         # write to files instead of screen/console
out.dir         <- "output/" # put generated files here

round.in.big.table <- 2

# Comment the next line if you want to keep old output files.
unlink(paste(out.dir, '*', sep=""))

# ----------   Don't touch anything past this line   ---------- #
# ---------- If you do, results might not replicate. ---------- #

# HELPERS

source('glmtools.R')


# SCRIPTS

cat("\014")

cat("\n\n LOADING ...\n")
source('1_load.R')       # imports and preparation/transformation of the data

cat("\n\n ATTRACTION ...\n")
source('2_attraction.R') # calculate attraction factors

cat("\n\n GLMM ...\n")
source('3_glmm.R')       # full GLMMs

#cat("\n\n GLMM with MCMC ...\n")
#source('4_glmm-mcmc.R')  # same GLMMs estimated with Stan (Bayesian/MCMC)

