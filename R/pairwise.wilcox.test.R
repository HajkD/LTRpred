pairwise.wilcox.test <- function(sim.matrix){
  
  pvals <- vector("numeric",ncol(sim.matrix)-2)
  for (i in 2:(ncol(sim.matrix)-1)){
    pvals[i - 1] <- stats::wilcox.test(x = sim.matrix[ , i],
                                       y = sim.matrix[ , i + 1],
                                       method = "two.sided", exact = FALSE)$p.value
    
    names(pvals)[i - 1] <- paste0(colnames(sim.matrix)[i],"<->",colnames(sim.matrix)[i + 1])
  }
  return (pvals)
}