require(plyr)


# Read data.
measure <- read.csv(file = "data/nan_sample.csv", head = TRUE, sep = '\t', quote = "")


# Backup original sentences.
sentences <- measure[, c("Leftcontext", "Match", "Rightcontext", "Originalsort", "Url")]


# Remove some columns. These are the ones we KEEP:
measure <- measure[,c("Measurecase",
          "Kindlemma",
          "Kindgender",
          "Kindcase",
          "Kindfreq",
          "Measurelemma",
          "Measureclass",
          "Measurefreq",
          "Minus1pos",
          "Genitives",
          "Badness",
          "Originalsort"
)]




# Rename some columns for better readability in paper.
colnames(measure)[which(colnames(measure) == "Minus1pos")] <- "Leftcontext"


# Filter data points that we cannot use.
measure <- measure[which(as.character(measure$Measurecase) %in% c("Nom","Acc","Dat")),]
measure <- measure[which(as.character(measure$Kindcase) %in% c("Acc", "Dat", "Nom", "Gen", "Obl", "Str", "Nomacc")),]
measure <- measure[-which(measure$Kindlemma %in% c("Sonnenschein","Silber","Obst","Strom")),]


# Simplify POS tags.
measure$Leftcontext <- as.character(measure$Leftcontext)
measure$Leftcontext <- ifelse(measure$Leftcontext %in% c("Art", "Pdat", "Pds", "Piat", "Pis", "Pposat"), "D", measure$Leftcontext)
measure$Leftcontext <- ifelse(measure$Leftcontext %in% c("Appr", "Apprart"), "P", measure$Leftcontext)
measure$Leftcontext <- ifelse(measure$Leftcontext == "Adja", "A", measure$Leftcontext)
measure$Leftcontext <- ifelse(measure$Leftcontext == "Card", "C", measure$Leftcontext)
measure$Leftcontext <- ifelse(measure$Leftcontext %in% c("D", "P", "A", "C"), measure$Leftcontext, "O")
measure$Leftcontext <- as.factor(measure$Leftcontext)

# Even simpler.
measure$Cardinal <- factor(ifelse(measure$Leftcontext == "C", "Yes", "No"), levels = c("Yes", "No"))

# Simplify "Measureclass".
measure$Measureclass <- factor(ifelse(as.character(measure$Measureclass) %in% c("Natural", "Trace", "Transport", "Object", "Unknown", "Unit", "Currency", "Layer"), "Rest", as.character(measure$Measureclass)))
measure$Measureclass <- factor(mapvalues(measure$Measureclass, from = c('Vessel', 'Length', 'Square', 'Capacity', 'Weight'), to = c('Container', rep('Physical', times = 4))))


# Remap abbreviated lemmas and some other lemma profliferation.
abbrs.from   <- c("mm", "ml", "kg", "qm", "DM", "cm", "l", "L", "km", "ha", "h", "g", "t", "mg", "€", "m", "kW", "T", "Kilo", "US-Dollar", "Leib", "Std.")
abbrs.to <- c("Millimeter", "Milliliter", "Kilogramm", "Quadratmeter", "Mark", "Zentimeter", "Liter", "Liter", "Kilometer", "Hektar", "Hektar", "Gramm", "Tonne", "Milligramm", "Euro", "Meter", "Kilowatt", "Tonne", "Kilogramm", "Dollar", "Laib", "Stunde")
measure$Measurelemma <- mapvalues(measure$Measurelemma, from = abbrs.from, to = abbrs.to)


# Remove uninformative cases (Fem kind with measure and kind both oblique).
measure <- measure[-which(measure$Kindgender == "Fem" & measure$Measurecase == "Dat"),]


# Now create response.
measure$Construction <- factor(ifelse(measure$Kindgender %in% c("Masc", "Neut"),
  ifelse(measure$Kindcase == "Gen", "PGCa", "NACa"),
  ifelse(measure$Kindcase == "Obl" & (measure$Measurecase == "Nom" | measure$Measurecase == "Acc"),"PGCa", "NACa")
  ), levels = c("NACa", "PGCa")
)


# Turn factors into factors and reorder them.
fax  <- c("Measurecase",
          "Kindlemma",
          "Kindgender",
          "Kindcase",
          "Measurelemma",
          "Measureclass",
          "Leftcontext")

for (fak in fax) {
  measure[,fak] <- factorify.to.zero(measure[,fak], measure$Genitive, T)
}


# Hm, different level ordering looks odd in plot, so make Leftcontext same order.
measure$Leftcontext  <- factor(measure$Leftcontext, levels = c("C", "D", "O", "A", "P"))
measure$Measureclass <- factor(measure$Measureclass, levels = c("Physical", "Container", "Rest", "Amount", "Portion" ))
measure$Kindgender   <- factor(measure$Kindgender, levels = c("Masc", "Neut", "Fem"))


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
  measure[,cntr] <- z.transform(measure[,cntr])
}

