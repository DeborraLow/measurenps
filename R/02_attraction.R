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

# Measure nousn have been pre-processed differently, can load in one table.
measurenouns_attracfreq <- read.delim("data/measurenouns_attracfreq.tsv", quote="")

# A cleanup.
colnames(ndn.kind) <- c("Freq","Lem")
colnames(nn.kind) <- c("Freq","Lem")
nn.kind <- nn.kind[-which(nn.kind$Lem=="rege"),]

# Get all lemmas for master table.
# Background: We want values for all nouns, even those which occur in
# only one of the NN/NDN tables. If not in one, count will just be 0.
cx.kind.lem <- c(as.character(nn.kind$Lem), as.character(ndn.kind$Lem))
cx.kind.lem <- unique(cx.kind.lem[order(cx.kind.lem)])

# Function to extract freqs for both constructions.
get_freq <- function(data, lem) {
  n1 <- data[which(data$Lem == lem),1]
  if (length(n1) == 0) n1 <- 0
  n1
}

# Generate full table.
cx.kind <- cbind(cx.kind.lem, lapply(cx.kind.lem, function(lem){get_freq(nn.kind,lem)}), lapply(cx.kind.lem, function(lem){get_freq(ndn.kind,lem)}))

# A rather simple "attraction" function, i.e., proportion.
attraction <- function(cs) {
  unlist(cs)[2] / (unlist(cs)[1] + unlist(cs)[2])
}



# Generate attraction quotients for full table.
cx.kind <- as.data.frame(cx.kind)
cx.kind$Pv <- unlist(apply(cx.kind[,2:3], 1, function(cs){attraction(cs)} ))

# Measure nouns are simpler because the tables were imported done & ready.
cx.measure <- measurenouns_attracfreq
cx.measure$Pv <- unlist(apply(cx.measure[,2:3], 1, function(cs){attraction(cs)} ))

# Give proper names to columns in final dfs.
colnames(cx.kind) <- c("Lemma", "Nn", "Ndn", "Pv")
colnames(cx.measure) <- c("Lemma", "Nn", "Ndn", "Pv")


# Pull attraction strength for each observation.

mn$Kindattraction <- unlist(lapply(mn$Kindlemma, function(x) { cx.kind[which(cx.kind$Lemma == x), "Pv"] } ))
mn$Measureattraction <- unlist(lapply(mn$Measurelemma, function(x) { cx.measure[which(as.character(cx.measure$Lemma) == as.character(x)), "Pv"] } ))
fem$Kindattraction <- unlist(lapply(fem$Kindlemma, function(x) { cx.kind[which(cx.kind$Lemma == x), "Pv"] } ))
fem$Measureattraction <- unlist(lapply(fem$Measurelemma, function(x) { cx.measure[which(as.character(cx.measure$Lemma) == as.character(x)), "Pv"] } ))

# Center
mn$Kindattraction <- center.simple(mn$Kindattraction)
mn$Measureattraction <- center.simple(mn$Measureattraction)
fem$Kindattraction <- center.simple(fem$Kindattraction)
fem$Measureattraction <- center.simple(fem$Measureattraction)
