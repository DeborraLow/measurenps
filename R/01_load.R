# 
#  LOAD DATA
#

require(gnumeric)
require(plyr)

linMap <- function(x, from, to)
  (x - min(x)) / max(x - min(x)) * (to - from) + from

m <- read.gnumeric.sheet(file = "data/nan_sample_KM.ods", head = TRUE, sheet.name = "Data")
pl <- read.gnumeric.sheet(file = "data/countnouns_sample_KM.ods", head = TRUE, sheet.name = "Data")

# Backup original sentences.

sentences <- m[, c("Leftcontext", "Match", "Rightcontext", "Originalsort", "Url")]

# Remove some columns. These are the ones we KEEP:

m <- m[,c("Measurecase",
          "Kindlemma",
          "Kindgender",
          "Kindcase",
          "Kindlength", 
          "Kindorigin",
          "Kindfinal",
          "Kindconsistency",
          "Kindedible",
          "Kindfreq",
          "Measurelemma",
          "Measuregender",
          "Measurenumber",
          "Measureabbreviated",
          "Measurelength",
          "Measureclass",
          "Measurefreq",
          "Minus1pos",
          "Minus2pos",
          "Matchlength",
          "Genitives",
          "Badness",
          "Originalsort"
)]

pl <- pl[,c("Measurecase",
          "Kindlemma",
          "Kindgender",
          "Kindcase",
          "Kindlength",
          # "Kindorigin",
          "Kindfinal",
          # "Kindconsistency",
          # "Kindedible",
          "Kindfreq",
          "Measurelemma",
          "Measuregender",
          "Measurenumber",
          "Measureabbreviated",
          "Measurelength",
          # "Measureclass",
          "Measurefreq",
          "Minus1pos",
          "Minus2pos",
          "Matchlength",
          "Genitives",
          "Badness",
          "Originalsort"
)]


# Filter data points that we cannot use.

m <- m[which(as.character(m$Measurecase) %in% c("Nom","Acc","Dat")),]
m <- m[which(as.character(m$Kindcase) %in% c("Acc", "Dat", "Nom", "Gen", "Obl", "Str", "Nomacc")),]
m <- m[-which(m$Kindlemma %in% c("Sonnenschein","Silber","Obst","Strom")),]
pl <- pl[which(as.character(pl$Measurecase) %in% c("Nom","Acc")),]
pl <- pl[which(as.character(pl$Kindcase) %in% c("Acc", "Dat", "Nom", "Gen", "Obl", "Str", "Nomacc")),]
pl <- pl[which(!as.character(pl$Measurelemma) %in% c("g", "l", "ml", "Kilo", "kg", "t", "Gramm", "Tonne", "Kilogramm")),]
pl$Measurelemma <- droplevels(pl$Measurelemma)


# Simplify POS tags.

m$Minus1pos <- as.character(m$Minus1pos)
m$Minus1pos <- ifelse(m$Minus1pos %in% c("Art", "Pdat", "Pds", "Piat", "Pis", "Pposat"), "D", m$Minus1pos)
m$Minus1pos <- ifelse(m$Minus1pos %in% c("Appr", "Apprart"), "P", m$Minus1pos)
m$Minus1pos <- ifelse(m$Minus1pos == "Adja", "A", m$Minus1pos)
m$Minus1pos <- ifelse(m$Minus1pos == "Card", "C", m$Minus1pos)
m$Minus1pos <- ifelse(m$Minus1pos %in% c("D", "P", "A", "C"), m$Minus1pos, "O")
m$Minus1pos <- as.factor(m$Minus1pos)

m$Minus2pos <- as.character(m$Minus2pos)
m$Minus2pos <- ifelse(m$Minus2pos %in% c("Art", "Pdat", "Pds", "Piat", "Pis", "Pposat"), "D", m$Minus2pos)
m$Minus2pos <- ifelse(m$Minus2pos %in% c("Appr", "Apprart"), "P", m$Minus2pos)
m$Minus2pos <- ifelse(m$Minus2pos == "Adja", "A", m$Minus2pos)
m$Minus2pos <- ifelse(m$Minus2pos == "Card", "C", m$Minus2pos)
m$Minus2pos <- ifelse(m$Minus2pos %in% c("D", "P", "A", "C"), m$Minus2pos, "O")
m$Minus2pos <- as.factor(m$Minus2pos)

pl$Minus1pos <- as.character(pl$Minus1pos)
pl$Minus1pos <- ifelse(pl$Minus1pos %in% c("Art", "Pdat", "Pds", "Piat", "Pis", "Pposat"), "D", pl$Minus1pos)
pl$Minus1pos <- ifelse(pl$Minus1pos %in% c("Appr", "Apprart"), "P", pl$Minus1pos)
pl$Minus1pos <- ifelse(pl$Minus1pos == "Adja", "A", pl$Minus1pos)
pl$Minus1pos <- ifelse(pl$Minus1pos == "Card", "C", pl$Minus1pos)
pl$Minus1pos <- ifelse(pl$Minus1pos %in% c("D", "P", "A", "C"), pl$Minus1pos, "O")
pl$Minus1pos <- as.factor(pl$Minus1pos)

pl$Minus2pos <- as.character(pl$Minus2pos)
pl$Minus2pos <- ifelse(pl$Minus2pos %in% c("Art", "Pdat", "Pds", "Piat", "Pis", "Pposat"), "D", pl$Minus2pos)
pl$Minus2pos <- ifelse(pl$Minus2pos %in% c("Appr", "Apprart"), "P", pl$Minus2pos)
pl$Minus2pos <- ifelse(pl$Minus2pos == "Adja", "A", pl$Minus2pos)
pl$Minus2pos <- ifelse(pl$Minus2pos == "Card", "C", pl$Minus2pos)
pl$Minus2pos <- ifelse(pl$Minus2pos %in% c("D", "P", "A", "C"), pl$Minus2pos, "O")
pl$Minus2pos <- as.factor(pl$Minus2pos)


# Create "Intentional" from "Origin".

m$Kindint <- factor(ifelse(m$Kindorigin %in% c("Inapplicable","Natural","Residue"), "Nonint", "Int"))


# Simplify "Measureclass".

m$Measureclass <- factor(ifelse(as.character(m$Measureclass) %in% c("Natural","Trace","Transport","Object","Unknown", "Unit"), "Rest", as.character(m$Measureclass)))



# Turn factors into factors and reorder them.

fax  <- c("Measurecase",
          "Kindlemma",
          "Kindgender",
          "Kindcase",
          "Kindfinal",
          "Kindconsistency",
          "Kindedible",
          "Kindorigin",
          "Kindint",
          "Measurelemma",
          "Measuregender", 
          "Measurenumber",
          "Measureabbreviated",
          "Measureclass",
          "Minus1pos",
          "Minus2pos")

no.pl.fax <- c(
  "Kindconsistency",
  "Kindedible",
  "Kindorigin",
  "Kindint",
  "Measureclass"
)


# Remap abbreviated lemmas and some other lemma profliferation.
abbrs.from   <- c("mm", "ml", "kg", "qm", "DM", "cm", "l", "L", "km", "ha", "h", "g", "t", "mg", "€", "m", "kW", "T", "Kilo", "US-Dollar", "Leib", "Std.")
abbrs.to <- c("Millimeter", "Milliliter", "Kilogramm", "Quadratmeter", "Mark", "Zentimeter", "Liter", "Liter", "Kilometer", "Hektar", "Hektar", "Gramm", "Tonne", "Milligramm", "Euro", "Meter", "Kilowatt", "Tonne", "Kilogramm", "Dollar", "Laib", "Stunde")
m$Measurelemma <- mapvalues(m$Measurelemma, from = abbrs.from, to = abbrs.to)
pl$Measurelemma <- mapvalues(pl$Measurelemma, from = abbrs.from, to = abbrs.to)



# Create binary factors.

m$Genitive <- factor(ifelse(m$Kindcase == "Gen", 1, 0), levels=c(0,1))
m$Casedrop <- factor(ifelse(m$Kindcase == "Obl" & (m$Measurecase == "Nom" | m$Measurecase == "Acc"), 1, 0), levels=c(0,1))
pl$Genitive <- factor(ifelse(pl$Kindcase == "Gen", 1, 0), levels=c(0,1))
pl$Casedrop <- factor(ifelse(pl$Kindcase == "Obl" & (pl$Measurecase == "Nom" | pl$Measurecase == "Acc"), 1, 0), levels=c(0,1))


# Split data by gender.

masc <- m[which(m$Kindgender=="Masc"),]
neut <- m[which(m$Kindgender=="Neut"),]
mn <- m[which(m$Kindgender=="Masc" | m$Kindgender=="Neut"),]
fem <- m[which(m$Kindgender=="Fem" & (m$Measurecase == "Nom" | m$Measurecase == "Acc")),]


# Push "Currency" TO "Rest" for Fem because there are too few feminine currencies.

fem$Measureclass <- factor(ifelse(as.character(fem$Measureclass) == "Currency", "Rest", as.character(fem$Measureclass)))

for (fak in fax) {
  mn[,fak] <- factorify.to.zero(mn[,fak], mn$Genitive, T)
}

for (fak in fax) {
  fem[,fak] <- factorify.to.zero(fem[,fak], fem$Casedrop, T)
}

for (fak in fax) {
  if (!fak %in% no.pl.fax) {
    print(fak)
    pl[,fak] <- factorify.to.zero(pl[,fak], pl$Casedrop, T)
  }
}

# Hm, different level ordering looks odd in plot, so make Minus1pos same order.
mn$Minus1pos <- factor(mn$Minus1pos, levels = c("C", "D", "O", "A", "P"))
fem$Minus1pos <- factor(fem$Minus1pos, levels = c("C", "D", "O", "A", "P"))


# Center numeric variables.

# This function centers mean to 0 and makes stdev:=1.
center.simple <- function(x) {
  x <- (x-mean(x))/sd(x)
  x
}

cntrs <- c("Kindfreq",
          "Measurelength",
          "Measurefreq",
          "Matchlength",
          "Genitives",
          "Badness")

for (cntr in cntrs) {
  mn[,cntr] <- center.simple(mn[,cntr])
  fem[,cntr] <- center.simple(fem[,cntr])
  pl[,cntr] <- center.simple(pl[,cntr])
}


