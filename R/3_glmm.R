require(lme4)
require(gridExtra)
require(MuMIn)
require(effects)
require(gridExtra)
require(car)


the.nAGQ               <- 0      # 0 for speed, 1 for precision and EXTREMELY slow computation.
ci.boot.nsim           <- 100   # Even 100 can be extremely slow! SET TO MIN. 1000 FOR PRODUCTION!!!
ci.boot.modelcomp.nsim <- 100   # This takes ages, even with 10. Set to 1000 for PRODUCTION!!!



measure.glmm <- glmer(Construction~1
                 
                 +(1|Measurelemma)
                 +(1|Kindlemma)
                 
                 +Badness                 
                 +Genitives
                 +Leftcontext
                 +Measurecase
                 +Measurenumber
                 
                 +Kindattraction
                 +Kindfreq
                 +Kindgender

                 +Measureattraction
                 +Measureclass
                 +Measurefreq

                 , data=measure, family=binomial(link=logit), na.action = na.fail, nAGQ=the.nAGQ,
                 control=glmerControl(optimizer="nloptwrap2", optCtrl=list(maxfun=2e5)))


# Standard diagnostics.
measure.glmm.wald <- Anova(measure.glmm, type="II")
measure.glmm.r2 <- r.squaredGLMM(measure.glmm)


# Predict and analyze with optimal cutoff.
all.cut           <- seq(-2,2,0.01)[which.max(unlist(lapply(seq(-2,2,0.01), function(x) {corr.prop(measure.glmm, measure$Construction, x)})))]
measure.glmm.pred <- ifelse(predict(measure.glmm) < all.cut, 0, 1)
measure.glmm.cmat <- table(measure.glmm.pred, measure$Construction, dnn=c("Pred","Obs"))
measure.glmm.corr <- sum(diag(measure.glmm.cmat))/sum(measure.glmm.cmat)
measure.glmm.base <- length(which(measure$Construction == "NACa"))/length(measure$Construction)
measure.glmm.pre  <- (measure.glmm.corr-measure.glmm.base)/(1-measure.glmm.base)


# Extract random effects.
measure.glmm.lemrand.order <- order(ranef(measure.glmm)$Kindlemma)
measure.glmm.lemrand       <- ranef(measure.glmm)$Kindlemma[measure.glmm.lemrand.order,]
measure.glmm.lemrand.names <- rownames(coef(measure.glmm)$Kindlemma)[measure.glmm.lemrand.order]



# PB test.
bootcomp.regs <- c("(1|Measurelemma)", "(1|Kindlemma)", 
                   "Badness", "Genitives", "Leftcontext", "Measurecase", "Measurenumber", 
                   "Kindattraction", "Kindfreq", "Kindgender", 
                   "Measureattraction", "Measureclass", "Measurefreq")
modelcomp <- lmer.modelcomparison(model = measure.glmm, regressors = bootcomp.regs,
                                     formula.target = "Construction~1", nsim = ci.boot.modelcomp.nsim,
                                     print.updated = T)



# Output.
if (save.persistent) sink(paste(out.dir, "glmm.txt", sep=""))
cat("\n\n\nMASCULINE\n\n")
print(summary(measure.glmm), correlation=F)
cat("\n\n")
print(measure.glmm.wald)
cat("\n\nModel comparison with LR and PB test\n")
print(modelcomp)
cat("\n\n")
print(measure.glmm.r2)
cat("\n\n")
cat("correct", measure.glmm.corr)
cat("\n\n")
cat("λ", measure.glmm.pre)
cat("\n\n")
if (save.persistent) sink()



# Get the bootstrapped CIs.
opts.ci.95 <- list(level = 0.95, method = "boot", boot.type = "perc", nsim = ci.boot.nsim, parallel="multicore", ncpus=8)
measure.ci.95 <- do.call(confint.merMod, c(opts.ci.95, list(object = measure.glmm, parm = names(fixef(measure.glmm)))))
measure.ci.95 <- measure.ci.95[nrow(measure.ci.95):2,]

opts.ci.90 <- list(level = 0.90, method = "boot", boot.type = "perc", nsim = ci.boot.nsim, parallel="multicore", ncpus=8)
measure.ci.90 <- do.call(confint.merMod, c(opts.ci.95, list(object = measure.glmm, parm = names(fixef(measure.glmm)))))
measure.ci.90 <- measure.ci.90[nrow(measure.ci.90):2,]

measure.fixeffs <- rev(fixef(measure.glmm)[2:length(fixef(measure.glmm))])

x.lower <- min(measure.ci.90[,1])*1.05
x.upper <- max(measure.ci.90[,2])*1.05

if (save.persistent) pdf(paste(out.dir, "fixeffs.pdf", sep=""))
dotchart(measure.fixeffs, pch=20,
         xlim = c(x.lower, x.upper),
         lcolor = "gray",
         cex = 1.2,
         main=paste("Fixed effects with bootstrapped\n 90% and 95% CIs (", ci.boot.nsim, " simulations)", sep=""))
lines(c(0,0), c(0,length(measure.ci.95)), col="gray")

for (i in 1:nrow(measure.ci.95)) {
  lines(measure.ci.90[i,c(1,2)], c(i,i), col="lightgray", lwd=5)
  points(measure.fixeffs[i], i, pch=18, cex=1.5, col="black")
  lines(measure.ci.95[i,c(1,2)], c(i,i), col="black", lwd=2)
}
if (save.persistent) dev.off()



# Effect plots.
fixeff.pl.cex <- 1.5
effs <- c("Badness", "Genitives", "Leftcontext", "Measurecase", "Measurenumber", 
          "Kindattraction", "Kindfreq", "Kindgender", 
          "Measureattraction", "Measureclass", "Measurefreq")

for (eff in effs) {
  if (save.persistent) pdf(paste(out.dir, eff, ".pdf", sep=""))
  p <- plot(effect(eff, measure.glmm),
            rug=F, 
            main = paste("EFfect plot: ", eff),
            ylab = "Probability of PGCa",
            colors = c("black", "darkblue"), cex = fixeff.pl.cex)
  print(p)
  if (save.persistent) dev.off()
}




# Selective ranef plots.
opts.dotchart <- list(pch=19, col="black", cex=1, xlab="Estimate of intercept")
n.select <- 30
main.dotchart.meas <- "Meas. rand. eff."
main.dotchart.kind <- "Kind rand. eff."

if (save.persistent) pdf(paste(out.dir, "raneffs_selection.pdf", sep=""))
par(mfrow=c(1,2))
do.call(ranef.plot, c(list(measure.glmm, measure, "Measurelemma", n.select, main = main.dotchart.meas), opts.dotchart))
do.call(ranef.plot, c(list(measure.glmm, measure, "Kindlemma", n.select, main = main.dotchart.kind), opts.dotchart))
par(mfrow=c(1,1))
if (save.persistent) dev.off()




# Make table.
measure.glmm.table <- data.frame(
  Coefficient = round(fixef(measure.glmm)[-1], round.in.big.table),
  Low         = rev(round(measure.ci.95[,1], round.in.big.table)),
  High        = rev(round(measure.ci.95[,2], round.in.big.table)),
  Not0        = rev(ifelse( (measure.ci.95[,1] < 0 & measure.ci.95[,2] < 0) | (measure.ci.95[,1] > 0 & measure.ci.95[,2] > 0) , "*", ""))
)

