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
setwd('~/Massangaben/')

# OPTIONS

set.seed(2398651)            # Do NOT change, or results will differ.
save.persistent <- T         # write to files instead of screen/console
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

if (save.persistent) sink(paste(out.dir, "00_wrapper.txt", ""), append=F)
cat("\nRoland Schäfer, Freie Universität Berlin\n")
cat("roland.schaefer@fu-berlin.de\n")
cat("http://orcid.org/0000-0003-3233-787\n")
cat("'Competing Constructions for German Measure Phrases'\n")
cat("Full results of corpus analysis as produced by R.\n")
cat("Run 'wrapper.R' for 100% reproducible results!")
cat("\n\nThis run BEGAN at ")
print(Sys.time())
if (save.persistent) sink()

# Steps 1 and 2 are required for the subsequent steps.

t.start <- Sys.time()

cat("\n\n LOADING ...\n")
source('01_load.R')       # imports and preparation/transformation of the data

t.loaded <- Sys.time()

cat("\n\n ATTRACTION ...\n")
source('02_attraction.R') # calculate attraction factors

t.attracted <- Sys.time()

cat("\n\n GLMM ...\n")
source('05_glmm.R')       # full GLMMs

t.glmm <- Sys.time()

cat("\n\n GLMM with MCMC ...\n")
source('06_glmm-mcmc.R')  # same GLMMs estimated with Stan (Bayesian/MCMC)

t.mcmc <- Sys.time()

if (save.persistent) sink(paste(out.dir, "00_wrapper.txt", ""), append=T)
cat("\n\nDONE.\n\n Time statistics:\n\n")

cat("Total:", t.mcmc - t.start)

cat("... Loading:", t.loaded - t.start, "\n")
cat("... Attraction:", t.attracted - t.loaded, "\n")
cat("... GLMMs:", t.glmm - t.attracted, "\n")
cat("... MCMC:", t.mcmc - t.glmm, "\n\n")

cat("\n\nThis run ENDED at ")
end.time <- Sys.time()
print(end.time)

if (save.persistent) sink()

# Save the state after all computations for debugging/ensuring reproducibility.
save(list = ls(), file=paste(out.dir, "workspace_", end.time, ".RData", sep=""))
