join.factors <- function(...) {
  unlist(list(...))
}

factorify <- function(f) {
  levelsize <- function(l) { s = length(which(f==l)) }
  f <- factor(f)
  o = levels(f)[order(unlist(lapply(levels(f), levelsize)), decreasing = T)]
  f <- factor(f, o)
}

# Sort factor to put most r(esponse)=0 in f(actor) first,
factorify.to.zero <- function(f, r, dec=F) {
  f <- factor(f)
  r <- factor(r)
  levs <- colnames(table(r, f))[order(table(r, f)[1,]/table(r, f)[2,],decreasing=dec)]
  f <- factor(f, levels=levs)
}


siglev <- function(p)
{
  if (p <= 0.1)
    if (p <= 0.05)
      if (p <= 0.01)
        if (p <= 0.001)
          sig <- "***"
  else
    sig <- "**"
  else
    sig <- "*"
  else
    sig <- "."
  else
    sig <- ""
}

siglev <- function(p)
{
  if (p <= 0.1)
    if (p <= 0.05)
      if (p <= 0.01)
        if (p <= 0.001)
          sig <- "***"
  else
    sig <- "**"
  else
    sig <- "*"
  else
    sig <- "."
  else
    sig <- ""
}

# This functions makes a plot of the COEFFICIENTS of significant factors.
plot.coef <- function(glm, alpha = 0.05, pr="Pr(>|z|)", cex.labels = 1, pos.labels = 0.4, ...)
{
  # Extract p values and coefs.
  p <- coef(summary(glm))[,pr]
  c <- coef(glm)
  
  # Separate intercept from other coefs.
  p.int <- as.numeric(p[1])
  c.int <- as.numeric(c[1])
  p <- p[-1]
  c <- c[-1]
  
  # Find the order in which to put coefficients.
  c <- c[which(p <= alpha)]
  p <- p[which(p <= alpha)]
  
  # Find the order according to coef.
  order <- order(c)
  c <- c[order]
  p <- p[order]
  names <- names(c)
  
  # Put intercept and coefs back together.
  names <- c("(Intercept)", names)
  vals <- c(round(c.int, 3), round(c, 3))
  
  c <- c(as.numeric(c.int), as.numeric(c))
  p <- c(as.numeric(p.int), as.numeric(p))
  vals <- paste(vals, lapply(p, siglev), sep=" ")
  names <- paste(names, vals, sep="\n")
  
  range <- max(c)-min(c)
  ylim <- c(min(c)-0.25*range, max(c)+0.25*range)
  bars <- barplot(c, ylim=ylim, col=c("white", rep("lightgray", length(c)-1)), ...)
  pos <- ifelse(c > 0, c+pos.labels*cex.labels, c-pos.labels*cex.labels)
  text(bars, pos, names, cex=cex.labels, srt=75)
}

# This functions makes a plot of the ODDS RATIOS of significant factors.
plot.or <- function(glm, alpha = 0.05, pr="Pr(>|z|)", cex.labels = 1, pos.labels = 3, ...)
{
  # Extract p values and coefs.
  p <- coef(summary(glm))[,pr]
  c <- exp(coef(glm))
  
  # Separate intercept from other coefs.
  p.int <- as.numeric(p[1])
  c.int <- as.numeric(c[1])
  p <- p[-1]
  c <- c[-1]
  
  # Extract the significant coefficients.
  c <- c[which(p <= alpha)]
  p <- p[which(p <= alpha)]
  
  # Find the order according to coef.
  order <- order(c)
  c <- c[order]
  p <- p[order]
  names <- names(c)
  
  # Put intercept and coefs back together.
  names <- c("(Intercept)", names)
  vals <- c(round(c.int, 3), round(c, 3))
  
  c <- c(as.numeric(c.int), as.numeric(c))
  p <- c(as.numeric(p.int), as.numeric(p))
  vals <- paste(vals, lapply(p, siglev), sep=" ")
  names <- paste(names, vals, sep="\n")
  
  range <- max(c)-min(c)
  ylim <- c(0, max(c)+0.25*range)
  bars <- barplot(c, ylim=ylim, col=c("white", rep("lightgray", length(c)-1)), ...)
  pos <- ifelse(c > 0, c+pos.labels*cex.labels, c-pos.labels*cex.labels)
  text(bars, pos, names, cex=cex.labels)
}


# Plots the intercept from model for the number most frequent levels in data.
ranef.plot <- function(model, data, effect, number, ...) {
  l.rle <- rle(as.character(data[order(as.character(data[,effect])),effect]))
  l.sel <- l.rle[["values"]][head(order(l.rle[["lengths"]], decreasing = T), number)]
  l.ranef.sel <- ranef(model)[[effect]][l.sel,]
  l.ranef.sel.order <- order(l.ranef.sel)
  dotchart(l.ranef.sel[l.ranef.sel.order], labels=l.sel[l.ranef.sel.order], ...)
}


lr.test <- function(glm, glm0)
{
  ll <- -2*logLik(glm)
  ll0 <- -2*logLik(glm0)
  lr <- ll0 - ll
  df <- glm$rank - glm0$rank
  p <- 1-pchisq(lr, df)
  list(lr=as.numeric(lr), df=as.numeric(df), p=as.numeric(p))
}

phi.glm <- function(glm)
{
  sum(resid(glm, type="pearson")^2 / df.residual(glm))
}
