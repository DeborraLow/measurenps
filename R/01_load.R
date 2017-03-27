# 
#  LOAD DATA
#

require(gnumeric)
require(plyr)

linMap <- function(x, from, to)
  (x - min(x)) / max(x - min(x)) * (to - from) + from

m <- read.gnumeric.sheet(file = "data/nan_sample_KM.ods", head = TRUE, sheet.name = "Data")

# Backup original sentences.

sentences <- m[, c("Leftcontext", "Match", "Rightcontext", "Originalsort", "Url")]

# Remove some columns. These are the ones we KEEP:

m <- m[,c("Measurecase",
          "Kindlemma",
          "Kindgender",
          "Kindcase",
          "Kindlength", 
          "Kindconsistency",
          "Kindedible",
          "Kindfreq",
          "Measurelemma",
          "Measurenumber",
          "Measureabbreviated",
          "Measureclass",
          "Measurefreq",
          "Minus1pos",
          "Genitives",
          "Badness",
          "Originalsort"
)]


# Filter data points that we cannot use.

m <- m[which(as.character(m$Measurecase) %in% c("Nom","Acc","Dat")),]
m <- m[which(as.character(m$Kindcase) %in% c("Acc", "Dat", "Nom", "Gen", "Obl", "Str", "Nomacc")),]
m <- m[-which(m$Kindlemma %in% c("Sonnenschein","Silber","Obst","Strom")),]


# Simplify POS tags.

m$Minus1pos <- as.character(m$Minus1pos)
m$Minus1pos <- ifelse(m$Minus1pos %in% c("Art", "Pdat", "Pds", "Piat", "Pis", "Pposat"), "D", m$Minus1pos)
m$Minus1pos <- ifelse(m$Minus1pos %in% c("Appr", "Apprart"), "P", m$Minus1pos)
m$Minus1pos <- ifelse(m$Minus1pos == "Adja", "A", m$Minus1pos)
m$Minus1pos <- ifelse(m$Minus1pos == "Card", "C", m$Minus1pos)
m$Minus1pos <- ifelse(m$Minus1pos %in% c("D", "P", "A", "C"), m$Minus1pos, "O")
m$Minus1pos <- as.factor(m$Minus1pos)


# Simplify "Measureclass".

m$Measureclass <- factor(ifelse(as.character(m$Measureclass) %in% c("Natural","Trace","Transport","Object","Unknown", "Unit", "Currency"), "Rest", as.character(m$Measureclass)))
m$Measureclass <- mapvalues(m$Measureclass, from = c('Vessel', 'Length', 'Square', 'Capacity', 'Weight'), to = c('Container', rep('Physical', times = 4)))



# Turn factors into factors and reorder them.

fax  <- c("Measurecase",
          "Kindlemma",
          "Kindgender",
          "Kindcase",
          "Kindconsistency",
          "Kindedible",
          "Measurelemma",
          "Measurenumber",
          "Measureabbreviated",
          "Measureclass",
          "Minus1pos")


# Remap abbreviated lemmas and some other lemma profliferation.
abbrs.from   <- c("mm", "ml", "kg", "qm", "DM", "cm", "l", "L", "km", "ha", "h", "g", "t", "mg", "€", "m", "kW", "T", "Kilo", "US-Dollar", "Leib", "Std.")
abbrs.to <- c("Millimeter", "Milliliter", "Kilogramm", "Quadratmeter", "Mark", "Zentimeter", "Liter", "Liter", "Kilometer", "Hektar", "Hektar", "Gramm", "Tonne", "Milligramm", "Euro", "Meter", "Kilowatt", "Tonne", "Kilogramm", "Dollar", "Laib", "Stunde")
m$Measurelemma <- mapvalues(m$Measurelemma, from = abbrs.from, to = abbrs.to)


# Create binary factors.

m$Genitive <- factor(ifelse(m$Kindcase == "Gen", 1, 0), levels=c(0,1))
m$Casedrop <- factor(ifelse(m$Kindcase == "Obl" & (m$Measurecase == "Nom" | m$Measurecase == "Acc"), 1, 0), levels=c(0,1))


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


# Hm, different level ordering looks odd in plot, so make Minus1pos same order.
mn$Minus1pos <- factor(mn$Minus1pos, levels = c("C", "D", "O", "A", "P"))
#mn$Measureclass <- factor(mn$Measureclass, levels = c("Length", "Square", "Capacity", "Weight", "Rest", "Currency", "Vessel", "Layer", "Amount", "Portion" ))
mn$Measureclass <- factor(mn$Measureclass, levels = c("Physical", "Rest", "Currency", "Container", "Layer", "Amount", "Portion" ))
fem$Minus1pos <- factor(fem$Minus1pos, levels = c("C", "D", "O", "A", "P"))




# Center numeric variables.

# This function centers mean to 0 and makes stdev:=1.
center.simple <- function(x) {
  x <- (x-mean(x))/sd(x)
  x
}

cntrs <- c("Kindfreq",
          "Measurefreq",
          "Genitives",
          "Badness")

for (cntr in cntrs) {
  mn[,cntr] <- center.simple(mn[,cntr])
  fem[,cntr] <- center.simple(fem[,cntr])
}


