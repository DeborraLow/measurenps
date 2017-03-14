library(plyr)

mn.selectors <- c("Minus1pos", "Measurenumber", "Measurecase")
prop.mn <- 0.15

fem.selectors <- c("Minus1pos", "Measurenumber", "Measurecase")
prop.fem <- 0.25

sel.no <- 7
per.cmb <- 2

exp_sample <- function(df, model, selectors, do.tail = F, prop = 0.25, find = 5, per.combo = 1, sel.feature = "Genitive") {
  pred.order <- order(predict(model))

  if (do.tail) {
    pred.osort <- df[tail(pred.order, prop*nrow(df)), "Originalsort"]
  } else {
    pred.osort <- df[head(pred.order, prop*nrow(df)), "Originalsort"]
  }

  exp.sel <- df[which(df$Originalsort %in% pred.osort), ]          # The actual data (part of big data frame).
  exp.conf <- ddply(exp.sel, selectors, nrow)                      # Just the feat/val configs of top pred cases.
  exp.conf  <- exp.conf[order(exp.conf$V1, decreasing = T),]       # The top summary/counting of feat/val configs.
  cat("Found feature configurations: ", nrow(exp.conf), " (showing 5)\n")
  print(head(exp.conf, 5))

  seen.knouns <- character()
  seen.mnouns <- character()
  
  found <- 0
  i <- 0
  
  # For each feature permutation, find unseen kind/measure lemmas.
  while (found < find) {
    i <- i + 1
        
    cat("\nFeatures selected: ", paste(as.character(unlist(exp.conf[i, 1:ncol(exp.conf)-1])), collapse="_"), "\n")
    
    if (do.tail) { sel.feature.val <- 1 } else { sel.feature.val <- 0 }
    tmp.sel <- exp.sel[which(exp.sel[,sel.feature] == sel.feature.val & unname(apply(exp.sel[,selectors], 1, function(x) { paste(as.character(x), collapse="_") })) == paste(as.character(unlist(exp.conf[i, 1:ncol(exp.conf)-1])), collapse="_")),]
    
    tmp.top.kind <- ddply(tmp.sel, "Kindlemma", nrow)
    tmp.top.knouns <- tmp.top.kind[order(tmp.top.kind$V1, decreasing = T),]
    tmp.top.meas <- ddply(tmp.sel, "Measurelemma", nrow)
    tmp.top.mnouns <- tmp.top.meas[order(tmp.top.meas$V1, decreasing = T),]
    
    tmp.tmp.sel <-tmp.sel[which(!(exp.sel$Measurelemma %in% seen.mnouns) & !(exp.sel$Kindlemma %in% seen.knouns)),]
    tmp.tmp.sel <- tmp.tmp.sel[which(complete.cases(tmp.tmp.sel)),]
    tmp.tmp.combs <- ddply(tmp.tmp.sel, c("Measurelemma", "Kindlemma"), nrow)
    tmp.tmp.combs <- tmp.tmp.combs[order(tmp.tmp.combs$V1, decreasing = T),]
    
    for (j in length(tmp.tmp.combs):1) {
      if (!(tmp.tmp.combs[j,"Measurelemma"] %in% seen.mnouns) & !(tmp.tmp.combs[j, "Kindlemma"] %in% seen.knouns)) {
        cat("   Lemmas selected: ", as.character(tmp.tmp.combs[j, "Measurelemma"]), as.character(tmp.tmp.combs[j, "Kindlemma"]), "\n")
        
        # TODO Get example.
        tmp.idx <- which(tmp.sel$Kindlemma == as.character(tmp.tmp.combs[j, "Kindlemma"]) & tmp.sel$Measurelemma == as.character(tmp.tmp.combs[j, "Measurelemma"]))

        if (length(tmp.idx) >= per.combo) {
          tmp.idx <- sample(tmp.idx, per.combo)
        } else {
          cat("           OUT OF LEMMAS! REUSING...\n")
          tmp.idx <- 1:nrow(tmp.sel)
          #if (length(tmp.idx) > 0)
          tmp.idx <- sample(tmp.idx, per.combo)
        }
        
        tmp.os <- tmp.sel[tmp.idx, "Originalsort"]

        for (k in 1:per.combo) {
          cat("      EX_L: ", as.character(sentences[which(sentences$Originalsort == tmp.os[k]),  "Leftcontext"]), "\n")
          cat("        _M: ", as.character(sentences[which(sentences$Originalsort == tmp.os[k]),  "Match"]), "\n")
          cat("        _R: ", as.character(sentences[which(sentences$Originalsort == tmp.os[k]),  "Rightcontext"]), "\n")          
        }
        
        seen.knouns[length(seen.knouns)+1] <- as.character(tmp.tmp.combs[j, "Kindlemma"])
        seen.mnouns[length(seen.mnouns)+1] <- as.character(tmp.tmp.combs[j, "Measurelemma"])
        found <-  found + 1
        break
      }
    }
  }
}

if (save.persistent) sink(paste(out.dir, "10_experiment.txt", sep = "", collapse = ""))

cat("\n\n SAMPLING STIMULI FOR EXPERIMENT \n\n")

cat("\n\n*** MASC ***\n")
cat("\n*** top  ***\n\n")
exp_sample(mn, mn.glmm, mn.selectors, F, prop.mn, sel.no, per.cmb, "Genitive")
cat("\n*** bot  ***\n\n")
exp_sample(mn, mn.glmm, mn.selectors, T, prop.mn, sel.no, per.cmb, "Genitive")

cat("\n")

cat("\n\n*** FEM ****\n")
cat("\n*** top  ***\n\n")
exp_sample(fem, fem.glmm, fem.selectors, F, prop.fem, sel.no, per.cmb, "Casedrop")
cat("\n*** bot  ***\n\n")
exp_sample(fem, fem.glmm, fem.selectors, T, prop.fem, sel.no, per.cmb, "Casedrop")

if (save.persistent) sink()


