require(lme4)
require(effects)
require(MuMIn)
require(pbkrtest)
require(plyr)
require(boot)
require(lattice)

rm(list = ls())

set.seed(239865)

setwd('/Users/user/Workingcopies/measurenps/GitHub/R')

target          <- "tf.rt.resid" # Response variable.
ci.method       <- "boot"        # boot
modcomp.nsim    <- 1000          # 1000
boot.nsim       <- 1000          # 1000
iqr.factor      <- 2             # For removing outliers in final data set.
out.dir         <- "output/spr_"
save.persistent <- T
data.dir        <- "data/spr/"
assigs.fn       <- "data/assignments.csv"
stims.fn        <- "data/spr_stimuli.csv"

data.files.names <- list.files(data.dir)
data.files.names <- data.files.names[grep('csv$', data.files.names)]
data.files       <- paste(data.dir, data.files.names, sep="")


load("output/corpus_stims_pred.RData")

assigs <- read.delim(assigs.fn, header = T, sep = ",", quote = "", stringsAsFactors = F)
stims  <- read.delim(stims.fn, header = T, sep = "\t", quote = "", stringsAsFactors = F)
stims$Modelprediction <- rep(c(inv.logit(stims.pred), rep(0, 33)), 2)

parts.rt.var <- NULL

# Get sentence lengths for stimuli (in order to get absolute indices into 
# reading time vectors later).
wlength <- function(x) length(unlist(strsplit(x, ' ', fixed = T)))
stims$Wordcount <- unlist(lapply(stims$Sentence, wlength))

rt.all <- data.frame()

# dafi <- data.files[1]
# Get data for each participant and add to master table.
for (dafi in data.files) {
  tmp <- read.delim(dafi, sep = ",", stringsAsFactors = F)
  
  # Get name, reading times, answers.
  part <- unname(unlist(tmp[1,"participant"]))
  rt   <- unname(unlist(tmp[which(tmp$wordRead.keys == "space"), "wordRead.rt"]))
  ans  <- unname(unlist(tmp[, "questAnswer.keys"]))
  
  print(part)
  
  # Get the stims for this part in part order and make right format (split words).
  part.stims.order   <- unname(unlist(assigs[which(assigs$Participant == part), 2:99]))
  part.stims         <- stims[part.stims.order, c("Sentence", "Target", "Kindgender", "Construction", "Modelpred", "Modelprediction", "Expectedcorr", "Adjtarget", "Wordcount", "Number", "Expno", "Item")]
  part.stims$woffset <- head(c(1, cumsum(part.stims$Wordcount)+1), -1)
  
  part.stims$Position = 1:nrow(part.stims)
  
  # Regress to get residual RTs.
  rt.resid.data           <- as.data.frame(cbind( nchar(unlist(strsplit(paste(part.stims$Sentence, collapse = " "), " ", fixed = T))), rt))
  rt.resid.data$rt.log    <- unlist(lapply(rt.resid.data$rt, log))
  colnames(rt.resid.data) <- c("wlength", "rt.raw", "rt.log")
  
  # Outlier correction.
  lo  <- quantile(rt.resid.data$rt.log)[2]-1.5*IQR(rt.resid.data$rt.log)
  hi  <- quantile(rt.resid.data$rt.log)[4]+1.5*IQR(rt.resid.data$rt.log)
  sel <- which(rt.resid.data$rt.log >= lo & rt.resid.data$rt.log <= hi)
  
  rt.resid.data.noutl <- rt.resid.data[sel,]
  
  # Residual model.
  rt.resid.model <- lm(rt.log~wlength, data = rt.resid.data.noutl)
  
  # Get basic data for the targets.
  part.stims$t1        <- apply(part.stims, 1, function(r) { ifelse(!r['Adjtarget'] == 'x', unlist(unname(strsplit(r["Sentence"], " ", fixed = T)))[as.numeric(r["Adjtarget"]) ], NA) } )
  part.stims$t2        <- apply(part.stims, 1, function(r) { ifelse(!r['Adjtarget'] == 'x', unlist(unname(strsplit(r["Sentence"], " ", fixed = T)))[as.numeric(r["Adjtarget"])+1 ], NA) } )
  part.stims$t1.length <- apply(part.stims, 1, function(r) { ifelse(!is.na(r["t1"]), nchar(r["t1"]), NA) } )
  part.stims$t2.length <- apply(part.stims, 1, function(r) { ifelse(!is.na(r["t2"]), nchar(r["t2"]), NA) } )
  part.stims$tf.length <- part.stims$t1.length + part.stims$t2.length

  # Get RTs.
  part.stims$t1.rt.raw <- apply(part.stims, 1, function(r) { ifelse(!is.na(r["t1"]), rt[as.numeric(r["woffset"])+as.numeric(r["Adjtarget"])-1], NA) } )
  part.stims$t2.rt.raw <- apply(part.stims, 1, function(r) { ifelse(!is.na(r["t2"]), rt[as.numeric(r["woffset"])+as.numeric(r["Adjtarget"])], NA) } )
  part.stims$tf.rt.raw <- part.stims$t1.rt.raw + part.stims$t2.rt.raw

  # Get log RTs.
  part.stims$t1.rt.log <- ifelse(!is.na(part.stims$t1.rt.raw), log(part.stims$t1.rt.raw), NA)
  part.stims$t2.rt.log <- ifelse(!is.na(part.stims$t2.rt.raw), log(part.stims$t2.rt.raw), NA)
  part.stims$tf.rt.log <- ifelse(!is.na(part.stims$tf.rt.raw), log(part.stims$tf.rt.raw), NA)

  # Get predicted RT by residual LM plus residual RT.
  part.stims$t1.rt.pred  <- apply(part.stims, 1, function(r) { predict(rt.resid.model, data.frame(wlength = as.numeric(r["t1.length"]))) } )
  part.stims$t1.rt.resid <- part.stims$t1.rt.log - part.stims$t1.rt.pred
  part.stims$t2.rt.pred  <- apply(part.stims, 1, function(r) { predict(rt.resid.model, data.frame(wlength = as.numeric(r["t2.length"]))) } )
  part.stims$t2.rt.resid <- part.stims$t2.rt.log - part.stims$t2.rt.pred
  part.stims$tf.rt.pred  <- apply(part.stims, 1, function(r) { predict(rt.resid.model, data.frame(wlength = as.numeric(r["tf.length"]))) } )
  part.stims$tf.rt.resid <- part.stims$tf.rt.log - part.stims$tf.rt.pred

  # Add participant label.
  part.stims$Participant <- part
  
  # Reduce to targets.
  part.stims <- part.stims[which(part.stims$Target==1),]
  
  if (is.null(rt.all)) {
    rt.all <- part.stims
  } else {
    rt.all <- rbind(rt.all, part.stims)
  }

  # DEBUG Plot relation per participant used for residualization.
  # plot(rt.resid.data.noutl$rt.log~rt.resid.data.noutl$wlength, pch=20, xlab="Wortlänge", ylab="log. Lesezeit", ylim = c(-3,2), main=paste(part))
  
  part.rt.avg      <- unname(apply(rt.resid.data.noutl, 1, function(r) { r["rt.raw"]/r["wlength"] }))
  part.rt.analysis <- data.frame(part = part, rt.mean = mean(part.rt.avg), rt.sd = sd(part.rt.avg))

  if (is.null(parts.rt.var)) {
    parts.rt.var <- part.rt.analysis
  } else {
    parts.rt.var <- rbind(parts.rt.var, part.rt.analysis)
  }
}



# Fix column types.
rt.all$Expectedcorr    <- as.factor(rt.all$Expectedcorr)
rt.all$Adjtarget       <- as.numeric(rt.all$Adjtarget)
rt.all$Modelpred       <- as.factor(rt.all$Modelpred)
rt.all$Construction    <- as.factor(rt.all$Construction)
rt.all$Kindgender      <- as.factor(rt.all$Kindgender)
rt.all$Modelprediction <- as.numeric(rt.all$Modelprediction)
rt.all$rt.tt.resid     <- as.numeric(rt.all$t1.rt.resid + rt.all$t1.rt.resid)


# Fix factor levels.
revalue(rt.all$Construction, c("NACa" = "NACadj", "PGCa" = "PGCadj"))


# ################
# Remove outliers.
# ################

before.oc      <- nrow(rt.all)
outlier.bounds <- as.numeric(c(quantile(unlist(rt.all[,target]))[3]-iqr.factor*IQR(unlist(rt.all[,target])),quantile(unlist(rt.all[,target]))[3]+iqr.factor*IQR(unlist(rt.all[,target]))))
rt.all         <- rt.all[which(rt.all[,target] >= outlier.bounds[1] & rt.all[,target] <= outlier.bounds[2]),]

if (save.persistent) sink(paste(out.dir, "results.txt", sep=""))
cat("\n\n Removed", before.oc-nrow(rt.all), "outliers between", outlier.bounds[1], "and", outlier.bounds[2], ".\n\n")
if (save.persistent) sink()

  
# ################################
# Make formula and estimate model.
# ################################

# Model with coprpus-based probs as input.
#formule.rsl <- paste(target, "~Construction*Modelprediction+Position+(1+Construction*Modelprediction|Participant)", sep = "")
formule <- paste(target, "~Construction*Modelprediction+Position+(1|Participant)+(1|Item)", sep = "")
formule.0 <- paste(target, "~Construction+Position+(1|Participant)+(1|Item)", sep = "")
#rt.model.rsl <- lmer(formule.rsl, data = rt.all, REML = T)
rt.model <- lmer(formule, data = rt.all, REML = T)
rt.model.0 <- lmer(formule.0, data = rt.all, REML = T)


# ######
# REPORT
# ######

if (save.persistent) sink(paste(out.dir, "results.txt", sep=""), append = T)

#cat("\n\n##### Model w/ random slope #####\n\n")
#print(summary(rt.model.rsl))
#cat("\n\n Model comparison \n")
#print(anova(rt.model, rt.model.0))
#cat("\n\n R-squared \n")
#cat("Full model\n")
#print(r.squaredGLMM(rt.model.rsl))

cat("\n\n##### Model w/o random slope #####\n\n")
print(summary(rt.model))
cat("\n\n Model comparison \n")
print(anova(rt.model, rt.model.0))
cat("\n\n R-squared \n")
cat("Full model\n")
print(r.squaredGLMM(rt.model))

cat("\n\n")
print(PBmodcomp(rt.model, rt.model.0), nsim = modcomp.nsim)

if (save.persistent) sink()


# #############
# Effect plots.
# #############

p <- plot(effect("Construction:Modelprediction", rt.model, KR = T), rug=F, colors = c("black", "darkblue"),
          ylab="Residual log. reading time", xlab = "Probability for PGCadj from corpus-based model",
          main=NULL)

if (save.persistent) {
  trellis.device(device = "pdf", file = paste0(out.dir, "effects.pdf"))
  trellis.par.set(list(axis.text = list(cex = 1.75)))
  trellis.par.set(list(par.ylab.text = list(cex = 1.75)))
  trellis.par.set(list(par.xlab.text = list(cex = 1.75)))
  trellis.par.set(list(par.zlab.text = list(cex = 1.75)))
  trellis.par.set(list(add.text = list(cex = 1.6)))
  trellis.par.set(list(par.sub.text = list(cex = 1.75)))
}
print(p)
if (save.persistent) dev.off()


# ###################
# Plot fixed effects.
# ###################

# Get bootstrap CIs 95%.
opts.ci.95 <- list(level = 0.95, method = ci.method, boot.type = "perc", nsim = boot.nsim)
rt.ci.95  <- do.call(confint.merMod, c(opts.ci.95, list(object = rt.model, parm = names(fixef(rt.model)))))
rt.ci.95  <- rt.ci.95[nrow(rt.ci.95):1,]  # Reverse order of CIs.

if (save.persistent) sink(paste(out.dir, "results.txt", sep=""), append = T)
cat("\n\n Bootstrap CIs\n")
print(rt.ci.95)
cat("\n\n")
if (save.persistent) sink()

fix.effs <- rev(fixef(rt.model)[1:length(fixef(rt.model))])

x.lower <- min(rt.ci.95[-nrow(rt.ci.95),1])*1.05
x.upper <- max(rt.ci.95[-nrow(rt.ci.95),2])*1.05

if (save.persistent) pdf(paste(out.dir, "fixefs.pdf", sep=""))
dotchart(fix.effs[-length(fix.effs)], pch=20,
         xlim = c(x.lower, x.upper),
         lcolor = "gray",
         cex = 0.8,
         main=paste("Fixed effects with bootstrapped\n 95% CIs (", boot.nsim, " simulations)", sep=""))
lines(c(0,0), c(0,length(rt.ci.95)), col="gray")

for (i in 1:(nrow(rt.ci.95)-1)) {
  points(fix.effs[i], i, pch=18, cex=1.5, col="black")
  lines(rt.ci.95[i,c(1,2)], c(i,i), col="black", lwd=2)
}
if (save.persistent) dev.off()






# #################
# Diagnostic plots.
# #################


if (save.persistent) pdf(paste(out.dir, "model_diagnostic.pdf", sep=""))
par(mfrow=c(3,2))
plot(residuals(rt.model) ~ unlist(rt.all[, target]), main="Model with probs", xlab="Reading time")
boxplot(residuals(rt.model) ~ rt.all$Construction, main="Model with probs", xlab="Case")
boxplot(residuals(rt.model) ~ rt.all$Modelprediction, main="Model with probs", xlab="Probs from corpus model")
boxplot(residuals(rt.model) ~ rt.all$Kindgender, main="Model with probs", xlab="Kind noun gender")
boxplot(residuals(rt.model) ~ rt.all$Position, main="Model with probs", xlab="Item position")
boxplot(residuals(rt.model) ~ rt.all$Participant, main="Model with probs", xlab="Participant")
par(mfrow=c(1,1))
if (save.persistent) dev.off()

if (save.persistent) pdf(paste(out.dir, "rt_participantvariance.pdf", sep=""))
lims <- c(min(parts.rt.var$rt.mean-parts.rt.var$rt.sd), max(parts.rt.var$rt.mean+parts.rt.var$rt.sd))
dotchart(parts.rt.var$rt.mean, pch="", xlim = lims,
         main = "Mean and sd character reading times\nper word for each participant"
         , labels = parts.rt.var$part
         #, labels = paste("Part", rev(1:nrow(parts.rt.var)), coll="")
         )
for (i in 1:nrow(parts.rt.var)) {
  lines(c(parts.rt.var[i, "rt.mean"]-parts.rt.var[i, "rt.sd"], parts.rt.var[i, "rt.mean"]+parts.rt.var[i, "rt.sd"]), c(i, i))
  points(parts.rt.var[i, "rt.mean"], i, pch = 18, cex = 1.5)
}
if (save.persistent) dev.off()



# #################
# GAMM
# #################

# To satisfy Divjak, Arppe & Baayen (2016) we also go fully GAMMy:

library(mgcv)

rt.gamm <- gamm(tf.rt.resid ~ s(Position) + s(Modelprediction, by = Construction),
  random = list( Participant = ~1 ),
  data = rt.all,
  family = gaussian,
  niterPQL = 1000
)

# ... with NO result. At least not different from LM.
# Notice that the relevant smoother (only significant term)
# is virtually linear.

if (save.persistent) sink(paste(out.dir, "results.txt", sep=""), append = T)
cat("\n\n GAMM, dedicated to Divjak, Arppe & Baayen (2016)\n")
print(summary(rt.gamm$gam))
if (save.persistent) sink()

if (save.persistent) pdf(paste(out.dir, "gamm.pdf", sep=""))
par(mfrow=c(2,2))
plot(rt.gamm$gam)
par(mfrow=c(1,1))
if (save.persistent) dev.off()


# Save workspace.
if (save.persistent) save(list = ls(), file=paste(out.dir, "workspace.RData", sep=""))
