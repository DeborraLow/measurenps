require(lme4)
require(gridExtra)
require(MuMIn)
require(effects)
require(gridExtra)
require(car)


the.nAGQ <- 0                 # 0 for speed, 1 for precision and EXTREMELY slow computation.
ci.boot.nsim <- 2            # Even 100 can be extremely slow! SET TO MIN. 1000 FOR PRODUCTION!!!
ci.boot.modelcomp.nsim <- 2  # This takes ages, even with 10. Set to 500 for PRODUCTION!!!



# Estimate model.
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

                 , data=m, family=binomial(link=logit), na.action = na.fail, nAGQ=the.nAGQ,
                 control=glmerControl(optimizer="nloptwrap2", optCtrl=list(maxfun=2e5)))


# Standard diagnostics.
measure.glmm.wald <- Anova(measure.glmm, type="II")
measure.glmm.r2 <- r.squaredGLMM(measure.glmm)


# Predict and analyze with optimal cutoff.
all.cut           <- seq(-2,2,0.01)[which.max(unlist(lapply(seq(-2,2,0.01), function(x) {corr.prop(measure.glmm, m$Construction, x)})))]
measure.glmm.pred <- ifelse(predict(measure.glmm) < all.cut, 0, 1)
measure.glmm.cmat <- table(measure.glmm.pred, m$Construction, dnn=c("Pred","Obs"))
measure.glmm.corr <- sum(diag(measure.glmm.cmat))/sum(measure.glmm.cmat)
measure.glmm.base <- length(which(m$Construction == "NACa"))/length(m$Construction)
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


# Effect plots.
fixeff.pl.cex <- 1.5
effs.all <- c("Genitives", "Leftcontext", "Measurenumber", "Badness", "Measurecase", "Kindattraction", "Measureattraction", "Measureclass",
             "Kindedible", "Kindfreq", "Measurefreq")

for (eff in effs.all) {
  print(eff)
  p <- plot(effect(eff, measure.glmm, KR = T), rug=F, main = paste("EFfect plot: ", eff), colors = c("black", "darkblue"), cex = fixeff.pl.cex)
  print(p)
}

