# This creates the quotient database for lemmas in the NN and NDN construction.
# It also creates some plots and a text file for illustration.

require(boot)

# How the plural lemma pair list was extracted from the filtered list in data frame pl:
# write.table(unique(sort(apply(pl[,c("Measurelemma", "Kindlemma")], 1, function(x) paste(x[1], x[2])))), "plural_pairs.txt", sep="\t", quote=F, col.names=F, row.names=F)

prop.func <- function(v, smooth = 1) {
  log( (v[2]+smooth) / (v[1]+smooth), 10)
}

ndn.kind <- read.delim("data/ndn_kind_freq.tsv", header=FALSE, quote="")
nn.kind <- read.delim("data/nn_kind_freq.tsv", header=FALSE, quote="")

# Pl has been pre-processed differently, can load in one table.
countnouns_attracfreq <- read.delim("data/countnouns_attracfreq.tsv", quote="")

# A cleanup.
nn.kind <- nn.kind[-which(nn.kind$Lem=="rege"),]
colnames(ndn.kind) <- c("Freq","Lem")
colnames(nn.kind) <- c("Freq","Lem")

# Get all lemmas for master table.
# Background: We want values for all nouns, even those which occur in
# only one of the NN/NDN tables. If not in one, count will just be 0.
cx.kind.lem <- c(as.character(nn.kind$Lem), as.character(ndn.kind$Lem))
cx.kind.lem <- unique(cx.kind.lem[order(cx.kind.lem)])

# Function to extract freqs for both constructions with "0 tolerance" (haha).
get_freq <- function(data, lem) {
  n1 <- data[which(data$Lem == lem),1]
  if (length(n1) == 0) n1 <- 0
  n1
}

# Generate full table.
cx.kind <- cbind(cx.kind.lem, lapply(cx.kind.lem, function(lem){get_freq(nn.kind,lem)}), lapply(cx.kind.lem, function(lem){get_freq(ndn.kind,lem)}))

# A rather simple "attraction" function, i.e., quotients.
attraction <- function(cs, smooth=1) {
  (unlist(cs)[2]+smooth) / (unlist(cs)[1]+smooth)
}

# Generate attraction quotients for full table.
cx.kind <- as.data.frame(cx.kind)
cx.kind$Pv <- unlist(apply(cx.kind[,2:3], 1, function(cs){attraction(cs)} ))
cx.kind$Lpv <- log(cx.kind$Pv, 10)

# PLURAL is simpler because the table was important done & ready.
cx.pl.kind <- countnouns_attracfreq
cx.pl.kind$Pv <- unlist(apply(cx.pl.kind[,2:3], 1, function(cs){attraction(cs)} ))
cx.pl.kind$Lpv <- log(cx.pl.kind$Pv, 10)

# Give proper names to columns in final dfs.
colnames(cx.kind) <- c("Lemma", "Nn", "Ndn", "Pv", "Lpv")
colnames(cx.pl.kind) <- c("Lemma", "Nn", "Ndn", "Pv", "Lpv")

# Order by attraction value.
cx.kind <- cx.kind[order (cx.kind$Pv),]
cx.pl.kind <- cx.pl.kind[order (cx.pl.kind$Pv),]

# Now, compare the "attraction" from neighboring constructions
# with the ratio Gen/CaseIdent ("prop") in the target constructions.

mn.props.kind <- apply(table(mn$Genitive,mn$Kindlemma), 2, prop.func )
fem.props.kind <- apply(table(fem$Casedrop,fem$Kindlemma), 2, prop.func )
pl.props.kind <- apply(table(pl$Casedrop,pl$Kindlemma), 2, prop.func )

# Now get the kind attration for the proportions from MN sample.
mn.assocs.kind <- unlist(lapply(names(mn.props.kind), function(x) { cx.kind[which(cx.kind$Lemma == x), "Lpv"] } ))
fem.assocs.kind <- unlist(lapply(names(fem.props.kind), function(x) { cx.kind[which(cx.kind$Lemma == x), "Lpv"] } ))
pl.assocs.kind <- unlist(lapply(names(pl.props.kind), function(x) { r <- cx.pl.kind[which(cx.pl.kind$Lemma == x), "Lpv"]; ifelse(is.null(r), 0, r) } ))

# Cleanup M/N.
mn.rmv <- which(names(mn.props.kind) %in% c("Regen","Wein","Blut","Vermögen"))
mn.nuller <- which(mn.props.kind=="0")
mn.outl <- c(mn.rmv,mn.nuller)
mn.props.kind.cl <- mn.props.kind[-mn.outl]
mn.assocs.kind.cl <- mn.assocs.kind[-mn.outl]

# Cleanup Fem.
fem.assocs.kind.cl <- fem.assocs.kind
fem.props.kind.cl <- fem.props.kind

# Cleanup Pl.
pl.assocs.kind.cl <- pl.assocs.kind
pl.props.kind.cl <- pl.props.kind

# Simple correlations and LMs.
mn.ass.corrtest <- cor.test(mn.props.kind.cl, mn.assocs.kind.cl)
mn.ass.lm <- lm(mn.props.kind.cl~mn.assocs.kind.cl)
fem.ass.corrtest <- cor.test(fem.props.kind.cl, fem.assocs.kind.cl)
fem.ass.lm <- lm(fem.props.kind.cl~fem.assocs.kind.cl)
pl.ass.corrtest <- cor.test(pl.props.kind.cl, pl.assocs.kind.cl)
pl.ass.lm <- lm(pl.props.kind.cl~pl.assocs.kind.cl)


# Bootstrap function wrapper for correlation.
boot.cor <- function(data, indices) {
  cor(data[indices,1], data[indices,2])
}

mn.cor.boot <- boot(data=cbind(mn.props.kind.cl, mn.assocs.kind.cl), statistic=boot.cor, R=100000, parallel="multicore", ncpus=8)
fem.cor.boot <- boot(data=cbind(fem.props.kind.cl, fem.assocs.kind.cl), statistic=boot.cor, R=100000, parallel="multicore", ncpus=8)
pl.cor.boot <- boot(data=cbind(pl.props.kind.cl, pl.assocs.kind.cl), statistic=boot.cor, R=100000, parallel="multicore", ncpus=8)


# OUTPUT

if (save.persistent) pdf(paste(out.dir, "02_mn_attraction_lm.pdf", sep=""))
plot(mn.props.kind.cl~mn.assocs.kind.cl)
abline(mn.ass.lm)
if (save.persistent) dev.off()

if (save.persistent) pdf(paste(out.dir, "02_fem_attraction_lm.pdf", sep=""))
plot(fem.props.kind.cl~fem.assocs.kind.cl)
abline(fem.ass.lm)
if (save.persistent) dev.off()

if (save.persistent) pdf(paste(out.dir, "02_pl_attraction_lm.pdf", sep=""))
plot(pl.props.kind.cl~pl.assocs.kind.cl)
abline(pl.ass.lm)
if (save.persistent) dev.off()

if (save.persistent) sink(paste(out.dir, "02_attraction.txt", sep=""))
cat("\n [ MASC / NEUT ATTRACTION EFFECT ]\n")
cat("\n [ Correlation ]\n\n")
print(mn.ass.corrtest)
cat("\n [ ... bootstrap CI ]\n\n")
print(mn.cor.boot)
print(boot.ci(mn.cor.boot, conf=0.95, type = "bca"))
cat("\n [ FEM ATTRACTION EFFECT ]\n")
cat("\n [ Correlation ]\n\n")
print(fem.ass.corrtest)
cat("\n [ ... bootstrap CI ]\n\n")
print(fem.cor.boot)
print(boot.ci(fem.cor.boot, conf=0.95, type = "bca"))
cat("\n [ PL ATTRACTION EFFECT ]\n")
cat("\n [ Correlation ]\n\n")
print(pl.ass.corrtest)
cat("\n [ ... bootstrap CI ]\n\n")
print(pl.cor.boot)
#print(boot.ci(pl.cor.boot, conf=0.95, type = "bca"))
if (save.persistent) sink()

if (save.persistent) {
  pdf(paste(out.dir, "02_mn_attraction_bootstrapcorr.pdf", sep=""))
  plot(mn.cor.boot)
  dev.off()
  pdf(paste(out.dir, "02_fem_attraction_bootstrapcorr.pdf", sep=""))
  plot(fem.cor.boot)
  dev.off()
#  pdf(paste(out.dir, "02_pl_attraction_bootstrapcorr.pdf", sep=""))
#  plot(pl.cor.boot)
#  dev.off()
}

# Pull attraction strength for each observation.

mn$Attraction <- unlist(lapply(mn$Kindlemma, function(x) { cx.kind[which(cx.kind$Lemma == x), "Lpv"] } ))
fem$Attraction <- unlist(lapply(fem$Kindlemma, function(x) { cx.kind[which(cx.kind$Lemma == x), "Lpv"] } ))
pl$Attraction <- unlist(lapply(as.character(pl$Kindlemma), function(x) { r <- cx.pl.kind[which(as.character(cx.pl.kind$Lemma) == x), "Lpv"]; ifelse(is.null(r), 0, r) } ))
pl$Attraction <- ifelse(is.na(pl$Attraction), 0, pl$Attraction)

# Center
mn$Attraction <- center.simple(mn$Attraction)
fem$Attraction <- center.simple(fem$Attraction)
pl$Attraction <- center.simple(pl$Attraction)
