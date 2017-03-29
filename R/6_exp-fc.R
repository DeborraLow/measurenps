require(lme4)
require(effects)
require(MuMIn)
require(pbkrtest)
require(boot)

rm(list = ls())

set.seed(239865)

setwd('/Users/user/Workingcopies/measurenps/GitHub/R')

ci.method        <- "boot" # boot
modcomp.nsim     <- 1000   # 1000
boot.nsim        <- 1000   # 1000
the.nAGQ         <- 0 
out.dir          <- "output/fc_"
save.persistent  <- T
data.dir         <- "data/fc/"
data.files.names <- list.files(data.dir)
data.files.names <- data.files.names[grep('csv$', data.files.names)]
data.files       <- paste(data.dir, data.files.names, sep="")

# Helpers.
source('glmtools.R')

load("output/corpus_stims_pred.RData")

results <- data.frame()

# Get "correctness" of responses.
for (dafi in data.files) {
  tmp <- read.delim(dafi, sep=",")
  tmp <- tmp[which(tmp$corrAns %in% c("l", "r")),]
  tmp$reaction.corr <- ifelse((tmp$reaction.keys == "f" & tmp$corrAns == "l") | (tmp$reaction.keys == "j" & tmp$corrAns == "r"), 1, 0)
  tmp <- tmp[order(tmp$trials.thisIndex),]
  results <- rbind(results, as.numeric(tmp$reaction.corr))
}

results <- t(results)

# This reconstructs the actual response from expected/not expected.
reconstruct.rating <- function(c) {
  !xor(rep(c(rep(0, 4), rep(1, 4)), 2), c)
}

# These are the predictions from the actual model.
responses.df <- data.frame(
    Modelpred     = rep(inv.logit(stims.pred), ncol(results))
  , Kindgender    = factor(rep(c(rep("MN", 8), rep("F", 8)), ncol(results)))
  , ResponseDrop  = factor(as.vector(apply(results, 2, reconstruct.rating)))
  , Item          = as.factor(rep(1:16, ncol(results)))
  , Participant   = as.factor(unlist(lapply(1:ncol(results), function(x) rep(x, 16))))
)


# Rename some variables to make paper consistent and readable.

colnames(responses.df)[which(colnames(responses.df) == "Modelpred")] <- "Modelprediction"
colnames(responses.df)[which(colnames(responses.df) == "ResponseDrop")] <- "Chosenconstruction"
responses.df$Chosenconstruction <- as.factor(ifelse(responses.df$Chosenconstruction == T, "PGCa", "NACa"))

# REPORT
if (save.persistent) sink(paste(out.dir, "results.txt", sep=""))

# Actual model.
cat("\n\nGLMM, predicting reactions from corpus model\n\n")
model.fc <- glmer(Chosenconstruction~Modelprediction+(1|Item)+(1|Participant),
                  family = binomial(link = "logit"), data = responses.df, nAGQ=the.nAGQ,
                  control=glmerControl(optimizer="nloptwrap2", optCtrl=list(maxfun=2e5)))
print(summary(model.fc))

# Model comparison.
cat("\n\nLR-Test and bootstrapped PB test on nested models w & w/o Modelpred\n\n")
model.fc.0 <- glmer(Chosenconstruction~(1|Item)+(1|Participant),
                  family = binomial(link = "logit"), data = responses.df, nAGQ=the.nAGQ,
                  control=glmerControl(optimizer="nloptwrap2", optCtrl=list(maxfun=2e5)))
print(PBmodcomp(model.fc, model.fc.0, nsim = modcomp.nsim))

cat("\n\n R-squared \n\n")
print(r.squaredGLMM(model.fc))

cat("\n\n 95% bootstrap CI")
opts.ci.95 <- list(level = 0.95, method = ci.method, boot.type = "perc", nsim = boot.nsim)
ci.95.fc <- do.call(confint.merMod, c(opts.ci.95, list(object = model.fc, parm = names(fixef(model.fc)))))
ci.95.fc <- ci.95.fc[nrow(ci.95.fc):1,]  # Reverse order of CIs and don't remove intercept. 

print(ci.95.fc)

if (save.persistent) sink()

# Effect plot for main interaction.
p <- plot(effect("Modelprediction", model.fc, KR = T), rug=F, colors = c("black", "darkblue"),
          main="Forced choice experiment",
          ylab="Probability that PGCa is chosen",
          xlab="Probability for PGCa from corpus-based model"
)
if (save.persistent) pdf(paste(out.dir, "effects.pdf", sep=""))
print(p)
if (save.persistent) dev.off()


# Plot coefs.
fix.effs <- rev(fixef(model.fc)[1:length(fixef(model.fc))])

x.lower <- min(ci.95.fc[,1])*1.05
x.upper <- max(ci.95.fc[,2])*1.05

if (save.persistent) pdf(paste(out.dir, "fixefs.pdf", sep=""))
dotchart(fix.effs, pch=20,
         xlim = c(x.lower, x.upper),
         lcolor = "gray",
         cex = 0.8,
         main=paste("Fixed effects with bootstrapped\n 95% CIs (", boot.nsim, " simulations)", sep=""))
lines(c(0,0), c(0,length(ci.95.fc)), col="gray")
for (i in 1:nrow(ci.95.fc)) {
  points(fix.effs[i], i, pch=18, cex=1.5, col="black")
  lines(ci.95.fc[i,c(1,2)], c(i,i), col="black", lwd=2)
}  
if (save.persistent) dev.off()


# Descriptive plot of responses model prediction.
if (save.persistent) pdf(paste(out.dir, "proportions.pdf", sep=""))
plot(responses.df$Chosenconstruction~responses.df$Modelprediction,
     main = "Forced choice: distribution of responses\nby (binned) predictions from corpus-based model",
     xlab = "Probability for PGCa from corpus-based model", ylab="Proportion of responses",
     col = c("gray30", "white"))
if (save.persistent) dev.off()


# Save workspace.
save(list = ls(), file=paste(out.dir, "workspace.RData", sep=""))

