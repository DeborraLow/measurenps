require(lme4)
require(gridExtra)
require(MuMIn)
require(effects)
require(gridExtra)


# ######### Some settings, conveniently at top of file. #########

the.nAGQ <- 0         # 0 for speed, 1 for precision and slow(!) computation.
ci.boot.nsim <- 1000  # Even 100 can be extremely slow! SET TO MIN. 1000 FOR PRODUCTION!!!
ci.boot.modelcomp.nsim <- 500 # This takes ages, even with 10. Set to 500 for PRODUCTION!!!

# ###############################################################


# Handling convergence warnings:
# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html

# A faster BOBYQA optimizer.
library(nloptr)
defaultControl <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1e-6,maxeval=1e5)
nloptwrap2 <- function(fn,par,lower,upper,control=list(),...) {
  for (n in names(defaultControl)) 
    if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
  res <- nloptr(x0=par,eval_f=fn,lb=lower,ub=upper,opts=control,...)
  with(res,list(par=solution,
                fval=objective,
                feval=iterations,
                conv=if (status>0) 0 else status,
                message=message))
}


# Function to get prop of correct predictions.
corr.prop <- function(model, observed, cutoff) {
  preds <- ifelse(predict(model) < cutoff, 0, 1)
  cmat <- table(preds, observed)
  sum(diag(cmat))/sum(cmat)
}


# Function to get LRT for (G)LMMMs with bootstrapped CIs. This takes AGES!
lmer.modelcomparison <- function(model, regressors, formula.target, ci = 0.95, nsim = 4, print.updated = F, print.diag = F) {
  require(pbkrtest)

  result = NULL
  
  for (i in 1:length(regressors)) {
    formula.parts <- c(formula.target, regressors[1:i-1], regressors[i+1:length(regressors)])
    formula.parts <- formula.parts[which(!is.na(formula.parts))]
    formula <- as.formula(paste(formula.parts, collapse="+"))
    
    model.0 <- update(model, formula = formula)
    pbmc <- PBmodcomp(model, model.0, nsim = nsim)
    
    if (print.diag) {
      cat("\n\n\n")
      print(regressors[i])
      print(pbmc)
    }
    
    res.tmp <- data.frame(Regressor = regressors[i], PB.p = pbmc$test$p.value[2],
                          LRT.p = pbmc$test$p.value[1], LR = unname(pbmc$LRTstat[1]),
                          Df = unname(pbmc$LRTstat[2]))
    if (is.null(result))
      result <- res.tmp
    else
      result <- rbind(result, res.tmp)
    
    if (print.updated) print(result)
  }
  result
}


# MASCULINE/NEUTER

# Note: High collinearity > 5 in Measureclass.
# But removing regressor does not change anything
# for the remaning factors substantially.
# Code:
# source("mer-utils.R")
# vif.mer(fem.glmm)

mn.glmm <- glmer(Genitive~1
                +(1|Measurelemma)
                +(1|Kindlemma)
                +Genitives
                +Minus1pos
                +Measureclass
                +Measureabbreviated
                +Measureattraction
                +Kindattraction
                +Measurenumber
                +Kindedible
                +Badness
                +Kindconsistency
                +Kindfreq
                +Measurecase
                +Measurefreq
              , data=mn, family=binomial(link=logit), na.action = na.fail, nAGQ=the.nAGQ,
              control=glmerControl(optimizer="nloptwrap2", optCtrl=list(maxfun=2e5)))


# Regressors to be used in getting LRT w/ bootstrapped CIs.

mn.bootcomp.regs <- c("(1|Measurelemma)", "(1|Kindlemma)", "Genitives", "Minus1pos",
                   "Measureclass", "Measureabbreviated", "Measureattraction",
                   "Kindattraction", "Measurenumber", "Kindedible",
                   "Badness", "Kindconsistency", "Kindfreq", "Measurecase", "Measurefreq")

# Get PB test.
mn.modelcomp <- lmer.modelcomparison(model = mn.glmm, regressors = mn.bootcomp.regs,
                     formula.target = "Genitive~1", nsim = ci.boot.modelcomp.nsim,
                     print.updated = T)


# Borin' ol' diagnostics.
mn.glmm.wald <- Anova(mn.glmm, type="II")
mn.glmm.r2 <- r.squaredGLMM(mn.glmm)

# Find best cutoff.
mn.cut <- seq(-2,2,0.01)[which.max(unlist(lapply(seq(-2,2,0.01), function(x) {corr.prop(mn.glmm, mn$Genitive, x)})))]
# Predict and analyze.
mn.glmm.pred <- ifelse(predict(mn.glmm) < mn.cut, 0, 1)
mn.glmm.cmat <- table(mn.glmm.pred, mn$Genitive, dnn=c("Pred","Obs"))
mn.glmm.corr <- sum(diag(mn.glmm.cmat))/sum(mn.glmm.cmat)
mn.glmm.base <- length(which(mn$Genitive == 0))/length(mn$Genitive)
mn.glmm.pre <- (mn.glmm.corr-mn.glmm.base)/(1-mn.glmm.base)
mn.glmm.lemrand.order <- order(ranef(mn.glmm)$Kindlemma)
mn.glmm.lemrand <- ranef(mn.glmm)$Kindlemma[mn.glmm.lemrand.order,]
mn.glmm.lemrand.names <- rownames(coef(mn.glmm)$Kindlemma)[mn.glmm.lemrand.order]

if (save.persistent) sink(paste(out.dir, "05_glmm.txt", sep=""))
cat("\n\n\nMASCULINE\n\n")
print(summary(mn.glmm), correlation=F)
cat("\n\n")
print(mn.glmm.wald)
cat("\n\nModel comparison with LR and PB test\n")
print(mn.modelcomp)
cat("\n\n")
print(mn.glmm.r2)
cat("\n\n")
cat("correct", mn.glmm.corr)
cat("\n\n")
cat("λ", mn.glmm.pre)
cat("\n\n")
if (save.persistent) sink()


# FEMININE

# Note: High collinearity > 5 in Kindconsistency.
# But removing regressor does not change anything
# for the remaning factors substantially.
# Code:
# source("mer-utils.R")
# vif.mer(fem.glmm)

fem.glmm <- glmer(Casedrop~1
                 +(1|Measurelemma)
                 +(1|Kindlemma)
                 +Measureclass
                 +Minus1pos
                 +Measureabbreviated
                 +Kindconsistency
                 +Measureattraction
                 +Genitives
                 +Kindedible
                 +Measurenumber
                 +Badness
                 +Measurefreq
                 +Measurecase
                 +Kindfreq
                 +Kindattraction
                 , data=fem, family=binomial(link=logit), nAGQ=the.nAGQ,
                 control=glmerControl(optimizer="nloptwrap2",optCtrl=list(maxfun=2e5)))

fem.bootcomp.regs <- c("(1|Measurelemma)", "(1|Kindlemma)", "Measureclass", "Minus1pos",
                       "Measureabbreviated", "Kindconsistency", "Measureattraction",
                       "Genitives", "Kindedible", "Measurenumber", "Badness", "Measurefreq",
                       "Measurecase", "Kindfreq", "Kindattraction")
# Get PB test.
fem.modelcomp <- lmer.modelcomparison(model = fem.glmm, regressors = fem.bootcomp.regs,
                     formula.target = "Casedrop~1", nsim = ci.boot.modelcomp.nsim,
                     print.updated = T, print.diag = T)

# Borin' ol' diagnostics.
fem.glmm.wald <- Anova(fem.glmm, type="II")
fem.glmm.r2 <- r.squaredGLMM(fem.glmm)
fem.cut <- seq(-2,2,0.01)[which.max(unlist(lapply(seq(-2,2,0.01), function(x) {corr.prop(fem.glmm, fem$Casedrop, x)})))]
fem.glmm.pred <- ifelse(predict(fem.glmm) < fem.cut, 0, 1)
fem.glmm.cmat <- table(fem.glmm.pred, fem$Casedrop)
fem.glmm.corr <- sum(diag(fem.glmm.cmat))/sum(fem.glmm.cmat)
fem.glmm.base <- length(which(fem$Casedrop == 0))/length(fem$Casedrop)
fem.glmm.pre <- (fem.glmm.corr-fem.glmm.base)/(1-fem.glmm.base)


if (save.persistent) sink(paste(out.dir, "05_glmm.txt", sep=""), append = T)
cat("\n\n\nFEMININE\n\n")
print(summary(fem.glmm, correlation=F))
cat("\n\n")
print(fem.glmm.wald)
cat("\n\nModel comparison with LR and PB test\n")
print(fem.modelcomp)
cat("\n\n")
print(fem.glmm.r2)
cat("\n\n")
cat("correct", fem.glmm.corr)
cat("\n\n")
cat("λ", fem.glmm.pre)
cat("\n\n")
if (save.persistent) sink()


# MORE OUTPUT

# Plot fixef and CI.

# Get the bootstrapped CIs.
opts.ci <- list(level = 0.95, method = "boot", boot.type = "perc", nsim = ci.boot.nsim, parallel="multicore", ncpus=8)
mn.ci <- do.call(confint.merMod, c(opts.ci, list(object = mn.glmm, parm = names(fixef(mn.glmm)))))
fem.ci <- do.call(confint.merMod, c(opts.ci, list(object = fem.glmm, parm = names(fixef(fem.glmm)))))

# Get Wald CIs.
opts.ci.wald <- list(level = 0.95, method = "Wald", boot.type = "perc", nsim = ci.boot.nsim, parallel="multicore", ncpus=8)
mn.ci.wald <- do.call(confint.merMod, c(opts.ci.wald, list(object = mn.glmm, parm = names(fixef(mn.glmm)))))
fem.ci.wald <- do.call(confint.merMod, c(opts.ci.wald, list(object = fem.glmm, parm = names(fixef(fem.glmm)))))


allfacs <- head(rev(c(names(fixef(mn.glmm)), names(fixef(fem.glmm))[!names(fixef(fem.glmm)) %in% names(fixef(mn.glmm))])), -1)
allfacs <- allfacs[rev(order(allfacs))]

col.fixef <- c("darkgreen", "lightgreen", "darkblue", "lightblue")
lwd.fixef <- c(3, 2)
pch.fixef <- c(25, 25, 24, 24)
cex.fixef <- 1.25

# This is a bit of a weird solution because
#  (a) it grew out of a simpler solution
#  (b) I should have encapsualted into a function but hadn't
for (zzz in 1:3) {

  if (zzz == 1) {
    localfacs <- allfacs[-grep('^(Kind|Measure)', allfacs)]
    localfacs <- c(localfacs, allfacs[grep('Measurecase|Measurenumber|Measureabbreviated', allfacs)])
    fnam  <- "04_glmm_fixef_firstlevel.pdf"
    subtitulum <- "First level regressors"
    posi <- "topright"
  }
  else if (zzz == 2) {
    localfacs <- allfacs[grep('^Measure', allfacs)]
    localfacs <- localfacs[which(!localfacs %in% allfacs[grep('Measurecase|Measurenumber|Measureabbreviated', allfacs)])]
    fnam  <- "04_glmm_fixef_secondlevel_measure.pdf"
    subtitulum <- "Second level regressors for measure lemma"
    posi <- "topright"
  }
  else {
    localfacs <- allfacs[grep('^Kind', allfacs)]
    fnam  <- "04_glmm_fixef_secondlevel_kind.pdf"
    subtitulum <- "Second level regressors for kind lemma"
    posi <- "topright"
  }

  localfacs <- sort(localfacs, decreasing = T)
  
  x.lower <- min(c(mn.ci[intersect(localfacs, rownames(mn.ci)),1], fem.ci[intersect(localfacs, rownames(fem.ci)),1]))*1.05
  x.upper <- max(c(mn.ci[intersect(localfacs, rownames(mn.ci)),2], fem.ci[intersect(localfacs, rownames(fem.ci)),2]))*1.05

  # Clamp for better readability of plot.
  x.lower <- ifelse(x.lower < -4, -4, x.lower)
  x.upper <- ifelse(x.upper > 4, 4, x.upper)
  
  if (save.persistent) pdf(paste(out.dir, fnam, sep=""))
  
  dotchart(rep(-10, length(localfacs)), xlim=c(x.lower, x.upper), lcolor = "gray", 
           labels = localfacs, cex = 0.7,
           main=paste("Fixed effects with bootstrapped CI:", subtitulum, sep="\n")
           )
  lines(c(0,0), c(0,length(localfacs)+1), col="gray", lty=1)
  
  for (i in 1:length(localfacs)) {
  
    # Draw F.
    if (localfacs[i] %in% names(fixef(fem.glmm))) {
      if ((fem.ci[localfacs[i],1] < 0 & fem.ci[localfacs[i],2] < 0) | (fem.ci[localfacs[i],1] > 0 & fem.ci[localfacs[i],2] > 0)) {
        lcol <- col.fixef[1]
        llty <- 1
        lpch <- pch.fixef[1]
        lcex <- cex.fixef
        llwd <- lwd.fixef[1]
      } else {
        lcol <- col.fixef[2]
        llty <- 1
        lpch <- pch.fixef[2]
        lcex <- cex.fixef
        llwd <- lwd.fixef[2]
      }
      if (fem.ci[localfacs[i],1] < x.lower) {
        x1 <- x.lower
        lines(c(x1-0.25, x1), c(i,i)+0.25, lty=3, col=lcol, cex = lcex, lwd = llwd)
      } else {
        x1 <- fem.ci[localfacs[i],1]
      }
      if (fem.ci[localfacs[i],2] > x.upper) {
        x2 <- x.upper
        lines(c(x2, x2+0.25), c(i,i)+0.25, lty=3, col=lcol, cex = lcex, lwd = llwd)
      } else {
        x2 <- fem.ci[localfacs[i],2]
      }
      lines(c(x1, x2), c(i,i)+0.25, lty=llty, lwd=llwd, col = lcol)
      points(fixef(fem.glmm)[localfacs[i]], i+0.25, col = lcol, bg = lcol, pch = lpch, cex = lcex)
    }  
    
    # Draw M/N.
    if (localfacs[i] %in% names(fixef(mn.glmm))) {
      if ((mn.ci[localfacs[i],1] < 0 & mn.ci[localfacs[i],2] < 0) | (mn.ci[localfacs[i],1] > 0 & mn.ci[localfacs[i],2] > 0)) {
        lcol <- col.fixef[3]
        llty <- 1
        lpch <- pch.fixef[3]
        lcex <- cex.fixef
        llwd <- lwd.fixef[1]
      } else {
        lcol <- col.fixef[4]
        llty <- 1
        lpch <- pch.fixef[4]
        lcex <- cex.fixef
        llwd <- lwd.fixef[2]
      }
      if (mn.ci[localfacs[i],1] < x.lower) {
        x1 <- x.lower
        lines(c(x1-0.25, x1), c(i,i)-0.25, lty=3, col=lcol, cex = lcex, lwd = llwd)
      } else {
        x1 <- mn.ci[localfacs[i],1]
      }
      if (mn.ci[localfacs[i],2] > x.upper) {
        x2 <- x.upper
        lines(c(x2, x2+0.25), c(i,i)-0.25, lty=3, col=lcol, cex = lcex, lwd = llwd)
      } else {
        x2 <- mn.ci[localfacs[i],2]
      }
      lines(c(x1, x2), c(i,i)-0.25, lty=llty, lwd=llwd, col = lcol)
      points(fixef(mn.glmm)[localfacs[i]], i-0.25, col = lcol, bg = lcol, pch = lpch, cex = lcex)
    }
  }
  
  legend(posi, bg = "white", legend = c("F (0 not in CI)", "F (0 in CI)", "M/N (0 not in CI)", "M/N (0 in CI)"), cex = 0.7,
         col = col.fixef, pt.bg = col.fixef, lty = 1, lwd = lwd.fixef, pch = pch.fixef, pt.cex = 0.65)
  
  if (save.persistent) dev.off()
}

# BIG TABLE (MCMC results to be added in next step)

allfacs <- rev(allfacs)
coefs.glmm.table <- data.frame(cbind(

  round(fixef(mn.glmm)[allfacs], round.in.big.table),
  round(mn.ci[,1][allfacs], round.in.big.table),
  round(mn.ci[,2][allfacs], round.in.big.table),
  round(mn.ci[,2][allfacs] - mn.ci[,1][allfacs], round.in.big.table),
  ifelse( (mn.ci[,1][allfacs] < 0 & mn.ci[,2][allfacs] < 0) | (mn.ci[,1][allfacs] > 0 & mn.ci[,2][allfacs] > 0) , "†", ""),
  ifelse( (mn.ci.wald[,1][allfacs] < 0 & mn.ci.wald[,2][allfacs] < 0) | (mn.ci.wald[,1][allfacs] > 0 & mn.ci.wald[,2][allfacs] > 0) , "†", ""),
  
  round(fixef(fem.glmm)[allfacs], round.in.big.table),
  round(fem.ci[,1][allfacs], round.in.big.table),
  round(fem.ci[,2][allfacs], round.in.big.table),
  round(fem.ci[,2][allfacs] - fem.ci[,1][allfacs], round.in.big.table),
  ifelse( (fem.ci[,1][allfacs] < 0 & fem.ci[,2][allfacs] < 0) | (fem.ci[,1][allfacs] > 0 & fem.ci[,2][allfacs] > 0) , "†", ""),
  ifelse( (fem.ci.wald[,1][allfacs] < 0 & fem.ci.wald[,2][allfacs] < 0) | (fem.ci.wald[,1][allfacs] > 0 & fem.ci.wald[,2][allfacs] > 0) , "†", "")
  
), row.names = allfacs)

colnames(coefs.glmm.table) <- c("MNCoef", "MNCIBootLo", "MNCIBootHi", "MNCIBootWidth", "MNCIBootEx0", "MNCIWald", "FCoef", "FCIBootLo", "FCIBootHi", "FCIBootWidth", "FCIBootEx0", "FCIWald")

if (save.persistent) sink(paste(out.dir, "05_glmm.txt", sep=""), append = T)
cat("\n\nTable comparing coefficient estimates between M/N and F (GLMM)\n\n")
print(coefs.glmm.table)
if (save.persistent) sink()


# Selective () ranef plots.
opts.dotchart <- list(pch=19, col="black", cex=0.8, xlab="Estimate of intercept")
n.select <- 30
main.dotchart.meas <- "Meas. rand. eff."
main.dotchart.kind <- "Kind rand. eff."

if (save.persistent) pdf(paste(out.dir, "04_glmm_raneff_mn.pdf", sep=""))
par(mfrow=c(1,1))
do.call(ranef.plot, c(list(mn.glmm, mn, "Kindlemma", n.select, main=paste(main.dotchart.kind, "(MN)")), opts.dotchart))
par(mfrow=c(1,1))
if (save.persistent) dev.off()

if (save.persistent) pdf(paste(out.dir, "04_glmm_raneff_fem.pdf", sep=""))
do.call(ranef.plot, c(list(fem.glmm, fem, "Measurelemma", n.select, main=paste(main.dotchart.meas, "(F)")), opts.dotchart))
if (save.persistent) dev.off()


# And full plots for personal use.

if (save.persistent) pdf(paste(out.dir, "04_glmm_raneff_mn_full.pdf", sep=""))
par(mfrow=c(1,1))
do.call(ranef.plot, c(list(mn.glmm, mn, "Kindlemma", length(unlist(ranef(mn.glmm)$Kindlemma)), main=paste(main.dotchart.kind, "(MN)")), opts.dotchart))
par(mfrow=c(1,1))
if (save.persistent) dev.off()

if (save.persistent) pdf(paste(out.dir, "04_glmm_raneff_fem_full.pdf", sep=""))
do.call(ranef.plot, c(list(fem.glmm, fem, "Measurelemma", length(unlist(ranef(fem.glmm)$Measurelemma)), main=paste(main.dotchart.meas, "(F)")), opts.dotchart))
if (save.persistent) dev.off()



# Effect plots.

fixeff.pl.cex <- 1.5

# M/N

effs.mn <- c("Genitives", "Minus1pos", "Measurenumber", "Badness", "Measurecase", "Kindattraction", "Measureattraction", "Measureclass",
             "Kindedible", "Kindfreq", "Measurefreq")

for (eff in effs.mn) {
  fn <- paste("output/04_glmm_fixeff_mn_", eff, ".pdf", sep="")
  if (save.persistent) pdf(fn)
  p <- plot(Effect(eff, mn.glmm, KR = T), rug=F, main = "Masculine/Neuter", colors = c("black", "darkblue"), cex = fixeff.pl.cex)
  print(p)
  if (save.persistent) dev.off()
}

# FEM

effs.fem <- c("Genitives", "Badness", "Measurecase", "Measureabbreviated", "Minus1pos", "Measurenumber",
              "Measureattraction", "Kindattraction", "Kindconsistency", "Kindedible", "Kindfreq", "Measurefreq")

for (eff in effs.fem) {
  print(eff)
  fn <- paste("output/04_glmm_fixeff_fem_", eff, ".pdf", sep="")
  if (save.persistent) pdf(fn)
  p <- plot(Effect(eff, fem.glmm, KR = T), rug=F, main = "Feminine", colors = c("black", "darkblue"), cex = fixeff.pl.cex,
            ylab = "Genitive")
  print(p)
  if (save.persistent) dev.off()
}



if (save.persistent) sink(paste(out.dir, "05_glmm.txt", sep=""), append = T)
source("mer-utils.R")
cat("\n\n MULTICOLLINEARITY\nM/N\n")
print(vif.mer(mn.glmm))
cat("\n\n F\n")
print(vif.mer(fem.glmm))
if (save.persistent) sink()

# Print raneffs,

if (save.persistent) sink(paste(out.dir, "05_glmm.txt", sep=""), append = T)
cat("\n\n Random effects\n")
cat("\n\n M/N Kindlemma \n")
cbind(
  rownames(ranef(mn.glmm)$Kindlemma)[order(ranef(mn.glmm)$Kindlemma[[1]])],
  unname(unlist(ranef(mn.glmm)$Kindlemma)[order(ranef(mn.glmm)$Kindlemma[[1]])])
)
cat("\n\n F Measurelemma \n")
cbind(
  rownames(ranef(fem.glmm)$Measurelemma)[order(ranef(fem.glmm)$Measurelemma[[1]])],
  unname(unlist(ranef(fem.glmm)$Measurelemma)[order(ranef(fem.glmm)$Measurelemma[[1]])])
)
if (save.persistent) sink()