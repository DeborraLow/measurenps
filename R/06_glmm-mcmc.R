require(rstanarm)
require(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# MASCULINE GLMM

chains <- 4
seed   <- 6976
iter   <- 1000 # Set to 1000 for production.

# MN

mn.glmm.mcmc <- stan_glmer(Genitive~1
                           +(1|Measurelemma)
                           +(1|Kindlemma)
                           +Genitives
                           +Minus1pos
                           #+Kindlength          # unstable in bootstrap
                           +Measureclass
                           +Matchlength
                           +Measureabbreviated
                           +Measureattraction
                           #+Measurelength       # virtual non-convergence
                           +Kindattraction
                           +Measurenumber
                           +Kindfinal
                           +Kindedible
                           +Badness
                           #+Kindconsistency     # virtual non-convergence
                           #+Kindorigin          # virtual non-convergence
                           +Kindfreq
                           +Measurecase
                           +Kindint
                           +Measurefreq
                           ,
                 data=mn, family=binomial(link=logit),
                 chains=chains, seed=seed, iter=iter,
                 prior = normal(0, 2.5), prior_intercept = normal(0, 10),
                 prior_ops = prior_options(prior_scale_for_dispersion = 5, min_prior_scale = 1e-12, scaled = TRUE),
                 prior_covariance = decov(regularization = 1, concentration = 1, shape = 1, scale = 1),
                 prior_PD = F
                 )

mn.baysum <- summary(mn.glmm.mcmc)


# FEM

fem.glmm.mcmc <- stan_glmer(Casedrop~1
                            +(1|Measurelemma)
                            +(1|Kindlemma)
                            #+Measureclass         # non-convergence
                            +Minus1pos
                            +Measureabbreviated
                            +Kindconsistency
                            +Matchlength
                            +Measureattraction
                            #+Kindfinal             # no hypothesis
                            +Genitives
                            +Kindedible
                            #+Kindorigin           # fixeff matrix rank deficient
                            #+Measurelength        # unstable in bootstrap
                            +Measuregender
                            +Measurenumber
                            +Kindint
                            +Badness
                            +Measurefreq
                            +Measurecase
                            +Kindfreq
                            +Kindattraction
                            #+Kindlength         # unstable in bootstrap
                  ,
                  data=fem, family=binomial(link=logit),
                  chains=chains, seed=seed, iter=iter,
                  prior = normal(0, 2.5), prior_intercept = normal(0, 10),
                  prior_ops = prior_options(prior_scale_for_dispersion = 5, min_prior_scale = 1e-12, scaled = TRUE),
                  prior_covariance = decov(regularization = 1, concentration = 1, shape = 1, scale = 1),
                  prior_PD = F
              )

fem.baysum <- summary(fem.glmm.mcmc)


# PL

pl.glmm.mcmc <- stan_glmer(Casedrop~1
                           +(1|Kindlemma)
                           #+(1|Measurelemma)    # model unidentifiable
                           #+Minus1pos           # model unidentifiable
                           +Attraction
                           #+Genitives           # model unidentifiable
                           #+Kindfinal           # model unidentifiable
                           #+Kindgender          # Hessian singular
                           #+Matchlength         # model unidentifiable
                           +Measurefreq
                           +Measurelength
                           #+Measurecase         # model unidentifiable
                           #+Kindfreq            # non-convergence
                           #+Badness             # model unidentifiable
                           #+Minus2pos           # model unidentifiable
                           #+Kindlength          # non-convergence
                           +Measuregender
                           #+Measurenumber       # makes no sense
                           #+Measureabbreviated  # only one level
                 ,
                 data=pl, family=binomial(link=logit),
                 chains=chains, seed=seed, iter=iter,
                 prior = normal(0, 2.5), prior_intercept = normal(0, 10),
                 prior_ops = prior_options(prior_scale_for_dispersion = 5, min_prior_scale = 1e-12, scaled = TRUE),
                 prior_covariance = decov(regularization = 1, concentration = 1, shape = 1, scale = 1),
                 prior_PD = F
               )

pl.baysum <- summary(pl.glmm.mcmc)



# OUTPUT

if (save.persistent) sink(paste(out.dir, "06_glmm-mcmc.txt", sep=""))
cat("\nBayesian estimation of GLMMs with MCMC\n\n")
print(mn.baysum)
print(fem.baysum)
print(pl.baysum)
if (save.persistent) sink()


# Add results to BIG TABLE.

coefs.mcmc.table <- data.frame(cbind(
  # M/N
  round(mn.baysum[,1][allfacs], round.in.big.table),
  round(mn.baysum[,4][allfacs], round.in.big.table),
  round(mn.baysum[,8][allfacs], round.in.big.table),
  round(mn.baysum[,8][allfacs] - mn.baysum[,4][allfacs], round.in.big.table),
  ifelse( (mn.baysum[,4][allfacs] < 0 & mn.baysum[,8][allfacs] < 0) | (mn.baysum[,4][allfacs] > 0 & mn.baysum[,8][allfacs] > 0) , "†", ""),
  # F
  round(fem.baysum[,1][allfacs], round.in.big.table),
  round(fem.baysum[,4][allfacs], round.in.big.table),
  round(fem.baysum[,8][allfacs], round.in.big.table),
  round(fem.baysum[,8][allfacs] - fem.baysum[,4][allfacs], round.in.big.table),
  ifelse( (fem.baysum[,4][allfacs] < 0 & fem.baysum[,8][allfacs] < 0) | (fem.baysum[,4][allfacs] > 0 & fem.baysum[,8][allfacs] > 0) , "†", "")
), row.names = allfacs)

colnames(coefs.mcmc.table) <- c("MNCoefMcmc", "MNCIMcmcLo", "MNCIMcmcHi", "MNCIMcmcWidth", "MNCIMcmcEx0", "FCoefMcmc", "FCIMcmcLo", "FCIMcmcHi", "FCIMcmcWidth", "FCIMcmcEx0")

# Create big table and order alphabetically.
big.table <- cbind(coefs.glmm.table[,1:6], coefs.mcmc.table[,1:5], coefs.glmm.table[,7:12], coefs.mcmc.table[,6:10])
big.table <- big.table[order(rownames(big.table)),]

if (save.persistent) sink(paste(out.dir, "06_glmm-mcmc.txt", sep=""), append = T)
cat("\n\nTable comparing coefficient estimates across relevant models\n\n")
print(big.table)
if (save.persistent) sink()