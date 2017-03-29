require(rstanarm)
require(rstan)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


measure.glmm.mcmc <- stan_glmer(Construction~1
                           
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
                           ,
                 data=measure, family=binomial(link=logit),
                 chains=chains, seed=seed, iter=iter,
                 prior = normal(0, 2.5), prior_intercept = normal(0, 10),
                 prior_ops = prior_options(prior_scale_for_dispersion = 5, min_prior_scale = 1e-12, scaled = TRUE),
                 prior_covariance = decov(regularization = 1, concentration = 1, shape = 1, scale = 1),
                 prior_PD = F
                 )

measure.baysum <- summary(measure.glmm.mcmc)


# OUTPUT

if (save.persistent) sink(paste(out.dir, "mcmc.txt", sep=""))
cat("\nBayesian estimation of GLMMs with MCMC\n\n")
print(measure.baysum)
if (save.persistent) sink()


# Make table.
measure.glmm.table <- data.frame(
  ML_Coef  = round(fixef(measure.glmm)[-1], precision),
  MC_Coef  = round(measure.baysum[,1][names(fixef(measure.glmm)[-1])], precision),
  
  ML_Low   = rev(round(measure.ci.95[,1], precision)),
  MC_Low   = round(measure.baysum[,4][names(fixef(measure.glmm)[-1])], precision),

  ML_High  = rev(round(measure.ci.95[,2], precision)),
  MC_High  = round(measure.baysum[,8][names(fixef(measure.glmm)[-1])], precision),

  ML_Not0  = rev(ifelse( (measure.ci.95[,1] < 0 & measure.ci.95[,2] < 0) | (measure.ci.95[,1] > 0 & measure.ci.95[,2] > 0) , "*", "")),
  MC_Not0  = ifelse( (measure.baysum[,4][names(fixef(measure.glmm)[-1])] < 0 & measure.baysum[,8][names(fixef(measure.glmm)[-1])] < 0) |
                       (measure.baysum[,4][names(fixef(measure.glmm)[-1])] > 0 & measure.baysum[,8][names(fixef(measure.glmm)[-1])] > 0) , "*", "")
)


# Print table.
if (save.persistent) sink(paste(out.dir, "mcmc.txt", sep=""), append = T)
cat("\n\nTable comparing coefficient estimates with ML and MCMC\n\n")
print(measure.glmm.table)
if (save.persistent) sink()


# Print new fixeff plot with MCMC.
if (save.persistent) pdf(paste(out.dir, "fixeffs_mle+mcmc.pdf", sep=""))
dotchart(rep(-100, length(measure.fixeffs)), pch=20,
         xlim = c(x.lower, x.upper),
         labels = names(measure.fixeffs),
         lcolor = "gray",
         cex = 1.2,
         main=paste("MLE and MCMC oefficient estimates\n with 95% confidence intervals", sep=""))
lines(c(0,0), c(0,length(measure.ci.95)), col="gray")

for (i in 1:nrow(measure.ci.95)) {
  points(measure.fixeffs[i], i+0.25, pch=25, cex=0.75, bg="black", col = "black")
  lines(measure.ci.95[i,c(1,2)], c(i+0.25, i+0.25), col="black", lwd=2)
  
  # MCMC.
  points(measure.glmm.table[names(measure.fixeffs[i]), "MC_Coef"], i-0.25, pch=24, cex=0.75, bg="gray", col = "gray")
  lines(measure.glmm.table[names(measure.fixeffs[i]), c("MC_Low", "MC_High")], c(i-0.25,i-0.25), col="gray", lwd=2)
}

legend("topright", legend = c("MLE", "MCMC"), cex = 1.2, col = c("black", "gray"), pt.bg = c("black", "gray"),
       lty = 1, lwd = 2, pch = c(25, 24), pt.cex = 0.75, bg = "white")
if (save.persistent) dev.off()