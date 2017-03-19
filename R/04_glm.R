require(fmsb)


###### MASCULINE ######

mn.glm <- glm(Genitive~1
              +Measurelemma
              #+Kindlemma       # prob 0 or 1 
              +Genitives
              +Minus1pos
              #+Kindlength       # unstable in bootstrap
              #+Measureclass     # prob 0 or 1
              +Matchlength
              +Measureabbreviated
              +Measureattraction
              +Measurelength
              +Kindattraction
              +Measurenumber
              +Kindedible
              +Badness
              +Kindconsistency
              +Kindorigin
              +Kindfreq
              #+Minus2pos       # prob 0 or 1 
              +Measurecase
              +Kindint
              +Kindgender
#              +Measurefreq
              , data=mn,  family=binomial(link=logit))

mn.glm.r2 <- NagelkerkeR2(mn.glm)$R2
mn.glm.pred <- ifelse(predict(mn.glm) < 0, 0, 1)
mn.glm.cmat <- table(mn.glm.pred, mn$Genitive)
mn.glm.corr <- sum(diag(mn.glm.cmat))/sum(mn.glm.cmat)
mn.glm.base <- length(which(mn$Genitive == 0))/length(mn$Genitive)
mn.glm.pre <- (mn.glm.corr-mn.glm.base)/(1-mn.glm.base)

if (save.persistent) sink(paste(out.dir, "04_glm.txt", sep=""))
cat("\n\n\nMASCULINE\n\n")
print(summary(mn.glm, correlation=F))
cat("\n\n")
print(mn.glm.r2)
cat("\n\n")
cat("correct", mn.glm.corr)
cat("\n\n")
cat("λ", mn.glm.pre)
cat("\n\n")
if (save.persistent) sink()



###### FEMININE ######


fem.glm <- glm(Casedrop~1
              +Measurelemma
              #+Kindlemma       # prob 0 or 1
              +Measureclass
              +Minus1pos
              +Measureabbreviated
              +Kindconsistency
              +Matchlength
              #+Measureattraction# prob 0 or 1
              +Genitives
              +Kindedible
              +Kindorigin
              +Measurelength
              +Measurenumber
              +Kindint
              +Badness
              +Minus1pos
              +Measurefreq
              +Measurecase
              +Kindfreq
              #+Kindattraction # prob 0 or 1
              #+Kindlength     # unstable in bootstrap
              , data=fem,  family=binomial(link=logit))
fem.glm.r2 <- NagelkerkeR2(fem.glm)$R2
fem.glm.pred <- ifelse(predict(fem.glm) < 0, 0, 1)
fem.glm.cmat <- table(fem.glm.pred, fem$Casedrop)
fem.glm.corr <- sum(diag(fem.glm.cmat))/sum(fem.glm.cmat)
fem.glm.base <- length(which(fem$Casedrop == 0))/length(fem$Casedrop)
fem.glm.pre <- (fem.glm.corr-fem.glm.base)/(1-fem.glm.base)

if (save.persistent) sink(paste(out.dir, "04_glm.txt", sep=""), append = T)
cat("\n\n\nFEMININE\n\n")
print(summary(fem.glm, correlation=F))
cat("\n\n")
print(fem.glm.r2)
cat("\n\n")
cat("correct", fem.glm.corr)
cat("\n\n")
cat("λ", fem.glm.pre)
cat("\n\n")
if (save.persistent) sink()


###### PLURAL ######


pl.glm <- glm(Casedrop~1
              +Kindlemma
              #+Measurelemma  # prob 0 or 1
              #+Minus1pos     # prob 0 or 1
              #+Attraction    # prob 0 or 1
              +Genitives
              #+Kindgender    # prob 0 or 1
              #+Matchlength   # prob 0 or 1
              #+Measurefreq   # prob 0 or 1
              #+Measurelength # prob 0 or 1
              +Measurecase
              #+Kindfreq      # prob 0 or 1
              +Badness
              #+Minus2pos     # prob 0 or 1
              +Kindlength
              #+Measurenumber # 1 level
              #+Measureabbreviated # 1 level
              , data=pl,  family=binomial(link=logit)
)
pl.glm.r2 <- NagelkerkeR2(pl.glm)$R2
pl.glm.pred <- ifelse(predict(pl.glm) < 0, 0, 1)
pl.glm.cmat <- table(pl.glm.pred, pl$Casedrop)
pl.glm.corr <- sum(diag(pl.glm.cmat))/sum(pl.glm.cmat)
pl.glm.base <- length(which(pl$Casedrop == 1))/length(pl$Casedrop)
pl.glm.pre <- (pl.glm.corr-pl.glm.base)/(1-pl.glm.base)

if (save.persistent) sink(paste(out.dir, "04_glm.txt", sep=""), append = T)
cat("\n\n\nPLURAL\n\n")
print(summary(pl.glm, correlation=F))
cat("\n\n")
print(pl.glm.r2)
cat("\n\n")
cat("correct", pl.glm.corr)
cat("\n\n")
cat("λ", pl.glm.pre)
cat("\n\n")
if (save.persistent) sink()
