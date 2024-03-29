
pairwise.z.test <- function(sim.matrix){
  
    sim.matrix <- as.data.frame(sim.matrix)
    pvals <- vector("numeric", ncol(sim.matrix) - 2)
    for (i in 2:(ncol(sim.matrix) - 1)) {
        pvals[i - 1] <- BSDA::z.test(
            x = sim.matrix[, i],
            y = sim.matrix[, i + 1],
            sigma.x = stats::sd(sim.matrix[, i]),
            sigma.y = stats::sd(sim.matrix[, i + 1]),
            alternative = "two.sided"
        )$p.value
        
        names(pvals)[i - 1] <-
            paste0(colnames(sim.matrix)[i], "<->", colnames(sim.matrix)[i + 1])
    }
    
    return (pvals)
  
}


  
  
