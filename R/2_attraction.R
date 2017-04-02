
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
cx.kind <- data.frame(
  Lemma = as.character(cx.kind.lem),
  Nn    = unlist(lapply(cx.kind.lem, function(lem){get_freq(nn.kind,lem)})),
  Ndn   = unlist(lapply(cx.kind.lem, function(lem){get_freq(ndn.kind,lem)}))
)


# A rather simple "attraction" function, i.e., log ratio.
attraction <- function(cs) {
  log((unlist(cs)[2] + 1) / (unlist(cs)[1] + 1))
}

collottraction <- function(cs, alls) {
  .cs <- as.numeric(cs)
  .alls <- as.numeric(alls)
  .mat<- matrix(c(.cs[2], .cs[1], .alls[2]-.cs[2], .alls[1]-.cs[1]), nrow = 2)
  .p <- fisher.test(.mat)$p.value
  log(.p, 10)
}



# Generate attraction quotients for full table.
kind.all.nn  <- sum(unlist((cx.kind$Nn)))
kind.all.ndn <- sum(unlist((cx.kind$Ndn)))
cx.kind <- as.data.frame(cx.kind)
cx.kind$Pv <- unlist(apply(cx.kind[,2:3], 1, function(cs){attraction(cs)} ))
cx.kind$Collo <- unlist(apply(cx.kind[,2:3], 1, function(cs){collottraction(cs, c(kind.all.nn, kind.all.ndn))} ))


# Measure nouns are simpler because the tables were imported done & ready.
cx.measure <- measurenouns_attracfreq
measure.all.nn  <- sum(unlist((cx.measure$Nn)))
measure.all.ndn <- sum(unlist((cx.measure$Ndn)))
cx.measure$Pv <- unlist(apply(cx.measure[,2:3], 1, function(cs){attraction(cs)} ))
cx.measure$Collo <- unlist(apply(cx.measure[,2:3], 1, function(cs){collottraction(cs, c(measure.all.nn, measure.all.ndn))} ))



# Pull attraction strength for each observation.
measure$Kindattraction <- unlist(lapply(measure$Kindlemma, function(x) { cx.kind[which(cx.kind$Lemma == as.character(x)), "Pv"] } ))
measure$Kindcollo <- unlist(lapply(measure$Kindlemma, function(x) { cx.kind[which(cx.kind$Lemma == as.character(x)), "Collo"] } ))
measure$Measureattraction <- unlist(lapply(measure$Measurelemma, function(x) { cx.measure[which(as.character(cx.measure$Lemma) == as.character(x)), "Pv"] } ))
measure$Measurecollo <- unlist(lapply(measure$Measurelemma, function(x) { cx.measure[which(as.character(cx.measure$Lemma) == as.character(x)), "Collo"] } ))

# Center
measure$Kindattraction    <- z.transform(measure$Kindattraction)
measure$Kindcollo         <- ifelse(measure$Kindcollo < -200, -200, measure$Kindcollo)
#measure$Kindcollo         <- z.transform(measure$Kindcollo)
measure$Measureattraction <- z.transform(measure$Measureattraction)
measure$Measurecollo      <- ifelse(measure$Measurecollo < -200, -200, measure$Measurecollo)
#measure$Measurecollo      <- z.transform(measure$Measurecollo)
