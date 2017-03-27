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