# Note: The pre-study was done AFTER the main study based on reviewer comments.
# This is why it is in a separate file and the coding style is slightly
# different.

require(pbkrtest)
require(lme4)
require(MuMIn)

rm(list = ls())

source('glmtools.R')

# z-transform numeric variables.
z.transform <- function(x) {
  x <- (x-mean(x))/sd(x)
  x
}

# Set this to the location of R script sources.
setwd('/Users/user/Workingcopies/measurenps/GitHub/R')

set.seed(953)
the.nAGQ <- 0
do.boot <- T
ci.boot.nsim <- 1000
ci.boot.modelcomp.nsim <- 1000
save.persistently <- T

mnps <- read.csv(file = "data/prestudy_concordance.csv", header = F, sep = "\t", quote = "",
                 colClasses = c(rep("factor", 5), "character", "factor", "factor", "numeric", rep("character", 3)),
                 col.names = c("Measurelemma", "Kindlemma", "Measureclass", "Construction", "Cardinal", "URL", "ID", "Badness",
                               "Genitives", "LC", "Match", "RC"))

kind <- read.csv(file = "data/prestudy_kinddata.csv", header = F, sep = "\t", quote = "",
                 colClasses = c("factor", "factor", "numeric"),
                 col.names = c("Kindlemma", "Kindgender", "Kindfreq"))

measure <- read.csv(file = "data/prestudy_measurefreq.csv", header = F, sep = "\t", quote = "",
                    colClasses = c("factor", "numeric"),
                    col.names = c("Measurelemma", "Measurefreq"))


# Transform badness and genitives as in main study.
mnps$Badness <- z.transform(sapply(as.character(mnps$Badness), utf8ToInt)-utf8ToInt('a'))


# Get the oldschool Genitive score from COW16 COReX genitive counts. See:
# https://github.com/rsling/cow/blob/master/src/de/old/cow-register-de
mnps$Genitives <- cut(mnps$Genitives, breaks = c(-1, 1, 2, 28, 41, 53, 65, 80, 102, 155, 1000))
levels(mnps$Genitives) <- as.character(seq(9, 0))
mnps$Genitives <- z.transform(as.numeric(levels(mnps$Genitives))[mnps$Genitives])

# Order factors in line with main study.
mnps$Measureclass <- factor(mnps$Measureclass, levels = c("Physical", "Container", "Rest", "Amount", "Portion"))
mnps$Cardinal <- factor(mnps$Cardinal, levels = c("Yes", "No"))

mnps <- merge(merge(mnps, kind, by = "Kindlemma"), measure, by = "Measurelemma")
mnps$Kindgender <- factor(mnps$Kindgender, levels = c("Masc", "Neut", "Fem"))
mnps$Measurefreq <- z.transform(mnps$Measurefreq)
mnps$Kindfreq <- z.transform(mnps$Kindfreq)

# Model.
mmodel <- glmer(Construction ~ Cardinal + Measureclass + Genitives + Badness + Measurefreq + Kindfreq  + Kindgender + (1|Measurelemma) + (1|Kindlemma),
                data = mnps, family=binomial(link=logit), na.action = na.fail, nAGQ=the.nAGQ,
                control=glmerControl(optimizer="nloptwrap2", optCtrl=list(maxfun=2e5)))
r2 <- r.squaredGLMM(mmodel)

if (do.boot) {
  # Bootstrap confidence.
  bootcomp.regs <- c("(1|Measurelemma)", "(1|Kindlemma)",
                     "Cardinal", "Measureclass", "Genitives", "Badness",
                     "Kindfreq", "Measurefreq", "Kindgender")
  modelcomp <- lmer.modelcomparison(model = mmodel, regressors = bootcomp.regs,
                                    formula.target = "Construction~1",
                                    nsim = ci.boot.modelcomp.nsim,
                                    print.updated = T)


  # Bootstrap model 'selection'.
  opts.ci.95 <- list(level = 0.95, method = "boot", boot.type = "perc", nsim = ci.boot.nsim, parallel="multicore", ncpus=2)
  ci95 <- do.call(confint.merMod, c(opts.ci.95, list(object = mmodel, parm = names(fixef(mmodel)))))
}

if (save.persistently) sink('output/prestudy.txt')
cat("=== MODEL SUMMARY ===\n\n")
print(summary(mmodel))
if (do.boot) {
  cat("\n\n=== BOOTSTRAP CONFIDENCE INTERVALS ===\n\n")
  print(ci95)
  cat("\n\n=== MODEL 'SELECTION ==='\n\n")
  print(modelcomp)
}
cat("\n\n=== R2 ===\n\n")
print(r2)
if (save.persistently) sink()




# New reduced model.
mmodel.redux <- glmer(Construction ~ Cardinal + Measureclass + Genitives + Badness + (1|Measurelemma) + (1|Kindlemma),
                data = mnps, family=binomial(link=logit), na.action = na.fail, nAGQ=the.nAGQ,
                control=glmerControl(optimizer="nloptwrap2", optCtrl=list(maxfun=2e5)))
r2.redux <- r.squaredGLMM(mmodel.redux)

if (do.boot) {
  # Bootstrap confidence.
  bootcomp.regs.redux <- c("(1|Measurelemma)", "(1|Kindlemma)",
                     "Cardinal", "Measureclass", "Genitives", "Badness")
  modelcomp.redux <- lmer.modelcomparison(model = mmodel.redux, regressors = bootcomp.regs.redux,
                                    formula.target = "Construction~1",
                                    nsim = ci.boot.modelcomp.nsim,
                                    print.updated = T)


# Bootstrap model 'selection'.
  ci95.redux <- do.call(confint.merMod, c(opts.ci.95, list(object = mmodel.redux, parm = names(fixef(mmodel)))))
}

if (save.persistently) sink('output/prestudy.txt', append = T)
cat("\n\n\n=== REDUCED MODEL ===\n\n")
cat("=== MODEL SUMMARY ===\n\n")
print(summary(mmodel.redux))
if (do.boot) {
  cat("\n\n=== BOOTSTRAP CONFIDENCE INTERVALS ===\n\n")
  print(ci95.redux)
  cat("\n\n=== MODEL 'SELECTION ==='\n\n")
  print(modelcomp.redux)
}
cat("\n\n=== R2 ===\n\n")
print(r2.redux)
if (save.persistently) sink()

