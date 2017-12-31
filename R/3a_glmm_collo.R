require(lme4)
require(gridExtra)
require(MuMIn)
require(effects)
require(gridExtra)
require(car)
require(pbkrtest)
require(lattice)

measure$Kindcollo <- z.transform(measure$Kindcollo)
measure$Measurecollo <- z.transform(measure$Measurecollo)

measure.glmm.collo <- glmer(Construction~1
                      
                      +(1|Measurelemma)
                      +(1|Kindlemma)
                      
                      +Badness                 
                      +Cardinal
                      +Genitives
                      +Measurecase

                      +Kindcollo
                      +Kindfreq
                      +Kindgender
                      
                      +Measurecollo
                      +Measureclass
                      +Measurefreq
                      
                      , data=measure, family=binomial(link=logit), na.action = na.fail, nAGQ=the.nAGQ,
                      control=glmerControl(optimizer="nloptwrap2", optCtrl=list(maxfun=2e5)))


# Standard diagnostics.
measure.glmm.collo.wald <- Anova(measure.glmm.collo, type="II")
measure.glmm.collo.r2 <- r.squaredGLMM(measure.glmm.collo)


# Predict and analyze with optimal cutoff.
all.cut           <- seq(-2,2,0.01)[which.max(unlist(lapply(seq(-2,2,0.01), function(x) {corr.prop(measure.glmm.collo, measure$Construction, x)})))]
measure.glmm.collo.pred <- ifelse(predict(measure.glmm.collo) < all.cut, 0, 1)
measure.glmm.collo.cmat <- table(measure.glmm.collo.pred, measure$Construction, dnn=c("Pred","Obs"))
measure.glmm.collo.corr <- sum(diag(measure.glmm.collo.cmat))/sum(measure.glmm.collo.cmat)
measure.glmm.collo.base <- length(which(measure$Construction == "NACa"))/length(measure$Construction)
measure.glmm.collo.pre  <- (measure.glmm.collo.corr-measure.glmm.collo.base)/(1-measure.glmm.collo.base)


# Extract random effects.
measure.glmm.collo.lemrand.order <- order(ranef(measure.glmm.collo)$Kindlemma)
measure.glmm.collo.lemrand       <- ranef(measure.glmm.collo)$Kindlemma[measure.glmm.collo.lemrand.order,]
measure.glmm.collo.lemrand.names <- rownames(coef(measure.glmm.collo)$Kindlemma)[measure.glmm.collo.lemrand.order]



# PB test.
set.seed(123)
model.0 <- update(measure.glmm.collo, formula = Construction~1+(1|Measurelemma)+(1|Kindlemma)+Badness+Cardinal+Genitives+
                    +Measurecase+Kindfreq+Kindgender+Measurecollo+Measureclass+Measurefreq)
pbmc.kindcollo <- PBmodcomp(measure.glmm.collo, model.0, nsim = ci.boot.modelcomp.nsim)
model.0 <- update(measure.glmm.collo, formula = Construction~1+(1|Measurelemma)+(1|Kindlemma)+Badness+Cardinal+Genitives+
                    +Measurecase+Kindfreq+Kindcollo+Kindgender+Measureclass+Measurefreq)
pbmc.measurecollo <- PBmodcomp(measure.glmm.collo, model.0, nsim = ci.boot.modelcomp.nsim)




# Output.
if (save.persistent) sink(paste(out.dir, "glmm_collo.txt", sep=""))
cat("\n\n\nGLMM with COLLO summary\n\n")
print(summary(measure.glmm.collo), correlation=F)
cat("\n\n")
print(measure.glmm.collo.wald)
cat("\n\nModel comparison with LR and PB test\n")
print(modelcomp)
print(round(modelcomp$PB.p, precision))
cat("\n\n")
print(round(measure.glmm.collo.r2, precision))
cat("\n\n")
cat("correct", round(measure.glmm.collo.corr, precision))
cat("\n\n")
cat("λ", round(measure.glmm.collo.pre, precision))
cat("\n\n")
if (save.persistent) sink()


