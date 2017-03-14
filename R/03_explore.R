#
# Explore by putting each potential regressor in simple GLM.
# ( Just to determine an order in which to try and add them to GLMM. )
#


# A function to do simple cross-validation.

crossval.glm <- function(model, K = 10) {
  data <- model$data
  response <- all.vars(model$formula)[1]
  n <- nrow(data)
  z <- n%/%K
  mod <- n%%K
  
  # Randomize order of data.
  idx <- sample(1:n, n, replace=F)
  data <- data[idx,]
  
  corrs <- c()
  
  # Iterate over the folds.
  for (i in 1:K) {
    
    if (i==K)
      idx.heldout <- seq((i-1)*z+1, i*z+mod, 1)
    else
      idx.heldout <- seq((i-1)*z+1, i*z, 1)
    
    data.train    <- data[-idx.heldout,]
    data.heldout  <- data[idx.heldout,]
    
    tryCatch({
      model.new <- glm(model$formula, family=model$family, data=data.train)
      preds <- ifelse(predict(model.new, newdata=data.heldout) < 0.5, 0, 1)
      confmat <- table(data.heldout[,response], preds)
      corr <- sum(diag(confmat))
      rate <- 1-corr/sum(confmat)
      corrs <- c(corrs, rate)
    },
    error = function(err) {
      cat("Hoppala! Convergence error...")
      # Yes, we just ignore errors! It's all about non-converging GLMs anyway.
    })
  }
  corr.armean <- 1/mean(1/corrs)
  corr.armean
}


# A function that performs a simple glm for each potential regressor.
# Should give a hint in which order to try and add regressors to complex model.

order.regressors <- function(df, target, ignores) {
  require(fmsb)
  require(boot)
  
  mat <- matrix(, nrow=length(names(df)), ncol=3)
  colnames(mat) <- c("R2", "p(LRT)", "Error(10CV)")
  rownames(mat) <- names(df)
  
  model.0 <- glm(df[,target]~1, family=binomial(link=logit))
  i = 1
  for (n in names(df)) {
    if (!(n %in% ignores)) {
      #cat(n, "\n") # To Debug models that dont work.
      
      formula <- paste(target, "~", n)
      model <- glm(as.formula(formula), family=binomial(link=logit), data=df)
      
      r2 <- NagelkerkeR2(model)
      lrt <- lr.test(model, model.0)
      delta <- crossval.glm(model)
      
      mat[i,] <- c(as.numeric(r2$R2), as.numeric(lrt$p), as.numeric(delta))
    } else {
      mat[i,] <- c(NA, NA, NA)
    }
    i = i + 1
  }
  mat <- mat[order(mat[,1], decreasing = T),]
  mat
}


# OUTPUT

if (save.persistent) sink(paste(out.dir, "03_explore.txt", sep=""))

ignores <- c("Kindcase", "Genitive", "Casedrop", "Originalsort")
mn.explore <- order.regressors(mn, "Genitive", ignores)
cat("MASCULINE\n")
print(rownames(mn.explore))

ignores <- c("Kindcase", "Genitive", "Casedrop", "Kindgender", "Originalsort")
fem.explore <- order.regressors(fem, "Casedrop", ignores)
cat("\n\nFEMININE\n")
print(rownames(fem.explore))

ignores <- c("Kindcase", "Genitive", "Casedrop", "Measurenumber", "Measureabbreviated", "Originalsort")
pl.explore <- order.regressors(pl, "Casedrop", ignores)
cat("\n\nPLURAL\n")
print(rownames(pl.explore))


# Create a full table.

# Calculate PRE (Kruskal's Lambda) for "explore" tables.
mn.baseline <- min(table(mn$Genitive))/sum(table(mn$Genitive))
fem.baseline <- min(table(fem$Casedrop))/sum(table(fem$Casedrop))
pl.baseline <- min(table(pl$Casedrop))/sum(table(pl$Casedrop))

mn.explore <- cbind(mn.explore, (mn.baseline-mn.explore[,3])/mn.baseline)
colnames(mn.explore) <- c(colnames(mn.explore)[1:3], "Lambda")
fem.explore <- cbind(fem.explore, (fem.baseline-fem.explore[,3])/fem.baseline)
colnames(fem.explore) <- c(colnames(fem.explore)[1:3], "Lambda")
pl.explore <- cbind(pl.explore, (pl.baseline-pl.explore[,3])/pl.baseline)
colnames(pl.explore) <- c(colnames(pl.explore)[1:3], "Lambda")

monofact.bigtable <- head(round(cbind(
  mn.explore[,c(1,2,4)], match(rownames(mn.explore), rownames(mn.explore)),
  fem.explore[rownames(mn.explore),c(1,2,4)], match(rownames(mn.explore), rownames(fem.explore)),
  pl.explore[match(rownames(mn.explore), rownames(pl.explore)),c(1,2,4)], match(rownames(mn.explore), rownames(pl.explore))
), 4), -4)

cat("\n\n Comparison of monofactorial effects for regressors.\n\n")
print(monofact.bigtable)

# Rank correlation between regressors ordered by monofactorial GLM.

cat("\n\nRank correlation between regressors ordered by monofactorial GLM.\n\n")

mnf.corr <- cor(order(rownames(mn.explore)[-1]), order(rownames(fem.explore)[-1]), method="sp")
cat('\n\nRank correlation (Spearman) MN/F : ', mnf.corr)

fem.explore.subs <- fem.explore[rownames(pl.explore),]
fem.explore.subs <- fem.explore.subs[order(fem.explore.subs[,1], decreasing = T),]
plf.corr <- cor(order(rownames(pl.explore)[-1]), order(rownames(fem.explore.subs)[-1]), method="sp")
cat('\n\nRank correlation (Spearman) PL/F : ', plf.corr)

mn.explore.subs <- mn.explore[rownames(pl.explore),]
mn.explore.subs <- mn.explore.subs[order(mn.explore.subs[,1], decreasing = T),]
plmn.corr <- cor(order(rownames(pl.explore)), order(rownames(mn.explore.subs)), method="sp")
cat('\n\nRank correlation (Spearman) PL/F : ', plf.corr)

if (save.persistent) sink()

if (save.persistent) pdf(paste(out.dir, "regressor_monofact_corr.pdf", sep=""))

plot(fem.explore[-1,1]~mn.explore[-1,1], pch=15, col="darkblue", xlim=c(0,0.45), ylim=c(0,0.45), cex=1.5,
     xlab="R²", ylab="R²", main="Correlation: monofactorial pre-analysis of regressors")
#abline(lm(fem.explore[-1,1]~mn.explore[-1,1]), col="darkblue", lty=1, lwd=2)
points(pl.explore[-1,1]~mn.explore.subs[-1,1], pch=16, col="orange", cex=1.5)
#abline(lm(pl.explore[-1,1]~mn.explore.subs[-1,1]), col="orange", lty=3, lwd=2)
points(fem.explore.subs[-1,1]~pl.explore[-1,1], pch=17, col="darkgreen", type="p", cex=1.5)
#abline(lm(fem.explore.subs[-1,1]~pl.explore[-1,1]), col="darkgreen", lty=2, lwd=2)
#abline(0,1, lwd=2, lty=2)
legend("bottomright", legend = c("Fem~Masc/Neut", "Pl~Masc/Neut", "Pl~Fem"), col=c("darkblue","orange","darkgreen"), pch=15:17)
#legend("topleft", legend = c("Fem~Masc/Neut", "Pl~Masc/Neut", "Pl~Fem"), col=c("darkblue","orange","darkgreen"), lty=1:3, pch=15:17)

if (save.persistent) dev.off()

# Check for multicolinearity.
source("highstat.r")

if (save.persistent) sink(paste(out.dir, "03_explore.txt", sep=""), append = T)

cat("\n\n Check for MULTICOLLINEARITY\n\n")

cat("\n == Masc/Neut ==\n")
corvif(mn[, c("Kindlength", "Kindfreq", "Measurelength", "Measurefreq", "Matchlength", "Genitives", "Badness", "Attraction")])

cat("\n == Fem ==\n")
corvif(fem[, c("Kindlength", "Kindfreq", "Measurelength", "Measurefreq", "Matchlength", "Genitives", "Badness", "Attraction")])

cat("\n == Pl ==\n")
corvif(pl[, c("Kindlength", "Kindfreq", "Measurelength", "Measurefreq", "Matchlength", "Genitives", "Badness", "Attraction")])

if (save.persistent) sink()