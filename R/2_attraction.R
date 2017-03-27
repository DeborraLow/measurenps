
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


# A rather simple "attraction" function, i.e., log ratio.
attraction <- function(cs) {
  log((unlist(cs)[2] + 1) / (unlist(cs)[1] + 1))
}


# Generate attraction quotients for full table.
cx.kind <- as.data.frame(cx.kind)
cx.kind$Pv <- unlist(apply(cx.kind[,2:3], 1, function(cs){attraction(cs)} ))


# Measure nouns are simpler because the tables were imported done & ready.
cx.measure <- measurenouns_attracfreq
cx.measure$Pv <- unlist(apply(cx.measure[,2:3], 1, function(cs){attraction(cs)} ))


# Give nice names to columns in final dfs.
colnames(cx.kind) <- c("Lemma", "Nn", "Ndn", "Pv")
colnames(cx.measure) <- c("Lemma", "Nn", "Ndn", "Pv")


# Pull attraction strength for each observation.
measure$Kindattraction <- unlist(lapply(measure$Kindlemma, function(x) { cx.kind[which(cx.kind$Lemma == x), "Pv"] } ))
measure$Measureattraction <- unlist(lapply(measure$Measurelemma, function(x) { cx.measure[which(as.character(cx.measure$Lemma) == as.character(x)), "Pv"] } ))


# Center
measure$Kindattraction    <- z.transform(measure$Kindattraction)
measure$Measureattraction <- z.transform(measure$Measureattraction)
