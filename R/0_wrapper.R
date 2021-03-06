
# Scripts for the analysis of the experiments and the pre-study need to be run manually.

rm(list = ls())


# Set this to the location of R script sources.
setwd('/Users/user/Workingcopies/measurenps/GitHub/R')


# Options.
set.seed(2398651)                          # Do NOT change, or results will differ.
save.persistent        <- T                # Write to files instead of screen/console.
out.dir                <- "output/"        # Put generated files here.
precision              <- 3                # Precision in result table.
the.nAGQ               <- 0                # 0 for speed, 1 for precision and eytremely slow computation.
ci.boot.nsim           <- 1000            # Even 100 can be extremely slow. Set to 1000 for production.
ci.boot.modelcomp.nsim <- 1000            # This takes ages, even with 10. Set to 1000 for for production.


# Delete old output.
if (save.persistent) unlink(paste(out.dir, '*', sep=""))
out.dir <- paste(out.dir, "corpus_", sep="")


# ----------   Don't touch anything past this line   ---------- #
# ---------- If you do, results might not replicate. ---------- #


# Helpers.
source('glmtools.R')

cat("\014")

cat("\n\n LOADING ...\n")
source('1_load.R')       # imports and preparation/transformation of the data

cat("\n\n ATTRACTION ...\n")
source('2_attraction.R') # calculate attraction factors

cat("\n\n GLMM ...\n")
source('3_glmm.R')       # full GLMMs

cat("\n\n GLMM with collexeme strength instead of attraction...\n")
source('3a_glmm_collo.R')       # full GLMMs will collonoscopy instead of attraction quotient.

cat("\n\n Prediction for experimental stimuli ...\n")
source('5_stimuli_predict.R')

# Save workspace.
if (save.persistent) save(list = ls(), file=paste(out.dir, "workspace.RData", sep=""))
