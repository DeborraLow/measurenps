require(gnumeric)
require(plyr)


# Read data.
m <- read.csv(file = "data/nan_sample.csv", head = TRUE, sep = '\t', quote = "")


# Backup original sentences.
sentences <- m[, c("Leftcontext", "Match", "Rightcontext", "Originalsort", "Url")]


# Remove some columns. These are the ones we KEEP:
m <- m[,c("Measurecase",
          "Kindlemma",
          "Kindgender",
          "Kindcase",
          "Kindfreq",
          "Measurelemma",
          "Measurenumber",
          "Measureclass",
          "Measurefreq",
          "Minus1pos",
          "Genitives",
          "Badness",
          "Originalsort"
)]




# Rename some columns for better readability in paper.
colnames(m)[which(colnames(m) == "Minus1pos")] <- "Leftcontext"


# Filter data points that we cannot use.
m <- m[which(as.character(m$Measurecase) %in% c("Nom","Acc","Dat")),]
m <- m[which(as.character(m$Kindcase) %in% c("Acc", "Dat", "Nom", "Gen", "Obl", "Str", "Nomacc")),]
m <- m[-which(m$Kindlemma %in% c("Sonnenschein","Silber","Obst","Strom")),]


# Simplify POS tags.
m$Leftcontext <- as.character(m$Leftcontext)
m$Leftcontext <- ifelse(m$Leftcontext %in% c("Art", "Pdat", "Pds", "Piat", "Pis", "Pposat"), "D", m$Leftcontext)
m$Leftcontext <- ifelse(m$Leftcontext %in% c("Appr", "Apprart"), "P", m$Leftcontext)
m$Leftcontext <- ifelse(m$Leftcontext == "Adja", "A", m$Leftcontext)
m$Leftcontext <- ifelse(m$Leftcontext == "Card", "C", m$Leftcontext)
m$Leftcontext <- ifelse(m$Leftcontext %in% c("D", "P", "A", "C"), m$Leftcontext, "O")
m$Leftcontext <- as.factor(m$Leftcontext)


# Simplify "Measureclass".
m$Measureclass <- factor(ifelse(as.character(m$Measureclass) %in% c("Natural","Trace","Transport","Object","Unknown", "Unit", "Currency"), "Rest", as.character(m$Measureclass)))
m$Measureclass <- mapvalues(m$Measureclass, from = c('Vessel', 'Length', 'Square', 'Capacity', 'Weight'), to = c('Container', rep('Physical', times = 4)))


# Remap abbreviated lemmas and some other lemma profliferation.
abbrs.from   <- c("mm", "ml", "kg", "qm", "DM", "cm", "l", "L", "km", "ha", "h", "g", "t", "mg", "€", "m", "kW", "T", "Kilo", "US-Dollar", "Leib", "Std.")
abbrs.to <- c("Millimeter", "Milliliter", "Kilogramm", "Quadratmeter", "Mark", "Zentimeter", "Liter", "Liter", "Kilometer", "Hektar", "Hektar", "Gramm", "Tonne", "Milligramm", "Euro", "Meter", "Kilowatt", "Tonne", "Kilogramm", "Dollar", "Laib", "Stunde")
m$Measurelemma <- mapvalues(m$Measurelemma, from = abbrs.from, to = abbrs.to)


# Remove uninformative cases (Fem kind with measure and kind both oblique).
m <- m[-which(m$Kindgender == "Fem" & m$Measurecase == "Dat"),]


# Now create response.
m$Construction <- factor(ifelse(m$Kindgender %in% c("Masc", "Neut"),
  ifelse(m$Kindcase == "Gen", "PGCa", "NACa"),
  ifelse(m$Kindcase == "Obl" & (m$Measurecase == "Nom" | m$Measurecase == "Acc"),"PGCa", "NACa")
  ), levels = c("NACa", "PGCa")
)


# Turn factors into factors and reorder them.
fax  <- c("Measurecase",
          "Kindlemma",
          "Kindgender",
          "Kindcase",
          "Measurelemma",
          "Measurenumber",
          "Measureclass",
          "Leftcontext")

for (fak in fax) {
  m[,fak] <- factorify.to.zero(m[,fak], m$Genitive, T)
}


# Hm, different level ordering looks odd in plot, so make Leftcontext same order.
m$Leftcontext <- factor(m$Leftcontext, levels = c("C", "D", "O", "A", "P"))
m$Measureclass <- factor(m$Measureclass, levels = c("Physical", "Rest", "Currency", "Container", "Layer", "Amount", "Portion" ))


# z-transform numeric variables.
z.transform <- function(x) {
  x <- (x-mean(x))/sd(x)
  x
}

cntrs <- c("Kindfreq",
          "Measurefreq",
          "Genitives",
          "Badness")

for (cntr in cntrs) {
  m[,cntr] <- z.transform(m[,cntr])
}

