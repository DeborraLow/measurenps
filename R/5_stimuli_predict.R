require(lme4)


stimuli <- read.csv("data/stimuli.csv", header=TRUE)


stims.fax  <- c(
  "Measurecase",
  "Kindlemma",
  "Measurelemma",
  "Measureclass",
  "Cardinal",
  "Kindgender"
  )


# Make factors factors, using levels from original model df.
for (fac in stims.fax) {
  stimuli[, fac]  <- factor(stimuli[, fac], levels = levels(measure[, fac]))
}


# Pull numeric values from original table because of standardization.
get.value <- function(df, index.column, index.value, value.column) {
  df[which(df[, index.column] == index.value)[1], value.column]
}

stimuli$Kindfreq          <- unlist(lapply(as.character(stimuli$Kindlemma), function(x) get.value(measure, "Kindlemma", x, "Kindfreq")))
stimuli$Measurefreq       <- unlist(lapply(as.character(stimuli$Measurelemma), function(x) get.value(measure, "Measurelemma", x, "Measurefreq")))
stimuli$Kindattraction    <- unlist(lapply(as.character(stimuli$Kindlemma), function(x) get.value(measure, "Kindlemma", x, "Kindattraction")))
stimuli$Measureattraction <- unlist(lapply(as.character(stimuli$Measurelemma), function(x) get.value(measure, "Measurelemma", x, "Measureattraction")))


# Set Badnesss to mean = 0 (b/o normalitzation).
stimuli$Badness <- 0
stimuli$Genitives <- 0


# Do the actual prediction.
stims.pred <- predict(measure.glmm, newdata = stimuli)

if (save.persistent) sink(paste(out.dir, "stimuli_predict.txt", sep=""))
print(stims.pred)
if (save.persistent) sink()

# Save to be loaded by experiment scripts.
save(stims.pred, file = paste(out.dir, "stims_pred.RData", sep=""))