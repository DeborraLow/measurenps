require(lme4)

stimuli <- read.csv("~/Massangaben/data/stimuli.csv", header=TRUE)

stims.fax  <- c("Measurecase",
          "Kindlemma",
          "Kindconsistency",
          "Kindedible",
          "Measurelemma",
          "Measurenumber",
          "Measureabbreviated",
          "Measureclass",
          "Minus1pos")


# Separate MN and F.

stimuli.mn <- stimuli[1:8,]
stimuli.fem <- stimuli[9:16,]


# Make factors factors, using levels from original model df.

for (fac in stims.fax) {
  stimuli.mn[, fac]  <- factor(stimuli.mn[, fac], levels = levels(mn[, fac]))
  stimuli.fem[, fac]  <- factor(stimuli.fem[, fac], levels = levels(fem[, fac]))
}


# Pull numeric values from original table because of standardization.

get.value <- function(df, index.column, index.value, value.column) {
  df[which(df[, index.column] == index.value)[1], value.column]
}

stimuli.mn$Kindfreq <- unlist(lapply(as.character(stimuli.mn$Kindlemma), function(x) get.value(mn, "Kindlemma", x, "Kindfreq")))
stimuli.mn$Measurefreq <- unlist(lapply(as.character(stimuli.mn$Measurelemma), function(x) get.value(mn, "Measurelemma", x, "Measurefreq")))
stimuli.mn$Kindattraction <- unlist(lapply(as.character(stimuli.mn$Kindlemma), function(x) get.value(mn, "Kindlemma", x, "Kindattraction")))
stimuli.mn$Measureattraction <- unlist(lapply(as.character(stimuli.mn$Measurelemma), function(x) get.value(mn, "Measurelemma", x, "Measureattraction")))

stimuli.fem$Kindfreq <- unlist(lapply(as.character(stimuli.fem$Kindlemma), function(x) get.value(fem, "Kindlemma", x, "Kindfreq")))
stimuli.fem$Measurefreq <- unlist(lapply(as.character(stimuli.fem$Measurelemma), function(x) get.value(fem, "Measurelemma", x, "Measurefreq")))
stimuli.fem$Kindattraction <- unlist(lapply(as.character(stimuli.fem$Kindlemma), function(x) get.value(fem, "Kindlemma", x, "Kindattraction")))
stimuli.fem$Measureattraction <- unlist(lapply(as.character(stimuli.fem$Measurelemma), function(x) get.value(fem, "Measurelemma", x, "Measureattraction")))


# Set Badnesss to mean = 0 (b/o normalitzation).
stimuli.mn$Badness <- 0
stimuli.mn$Genitives <- 0
stimuli.fem$Badness <- 0
stimuli.fem$Genitives <- 0


# Do the actual prediction.

if (save.persistent) sink(paste(out.dir, "11_Stimuli_predict.txt", sep=""))

print(predict(mn.glmm, newdata = stimuli.mn))
print(predict(fem.glmm, newdata = stimuli.fem))

if (save.persistent) sink()

# Check whether all variables were correctly transformed.
# for (s in 7:ncol(stimuli.mn)) {
#   coln <- colnames(stimuli.mn[s])
#   cat("\n\n\n\n", coln, "\n\n")
#   cat("Model\n")
#   print(summary(mn[, coln]))
#   cat("New\n")
#   print(summary(stimuli.mn[, coln]))  
# }