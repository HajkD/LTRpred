#' @title Plot LTR Similarity vs. predicted LTR count
#' @description The LTR Similarity vs. predicted LTR count visualization
#' allows to study the variability of predicted LTRs in different genomes.
#' @param sim.matrix Similarity matrix retruned by \code{\link{LTRpred.meta}}.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param main main text.
#' @param text.size size of the labels in ggplot2 notation.
#' @author Hajk-Georg Drost
#' @examples 
#' SimMatrix <- read.csv(system.file("SimMatrix.csv",package = "LTRpred"), sep = ";")
#' GenomeMatrix <- read.csv(system.file("GenomeMatrix.csv",package = "LTRpred"), sep = ";")
#' names(SimMatrix) <- c("organism",levels(cut(100,rev(seq(100,70,-2)),
#'                                  include.lowest = TRUE,right = TRUE)))
#' SimMatrix2 <- SimMatrix
#' SimMatrix2[ , 2:16] <- SimMatrix2[ , 2:16] / GenomeMatrix$genome.size
#' 
#' PlotSimCount(SimMatrix)
#' PlotSimCount(SimMatrix2)
#' @export

PlotSimCount <- function(sim.matrix, 
                         xlab            = "LTR % Similarity", 
                         ylab            = "Count", 
                         main            = "LTR % Similarity vs. Count",
                         text.size       = 18){
  
  pvals <- pairwise.z.test(sim.matrix)
  
  pvals.name <- rep("",ncol(sim.matrix)-1)
  pvals.name[which(pvals <= 0.05)] <- "*"
  pvals.name[which(pvals <= 0.005)] <- "**"
  pvals.name[which(pvals <= 0.0005)] <- "***"
  pvals.name[which(is.na(pvals))] <- ""
  
  colnames(sim.matrix)[1] <- "organism"
  reshaped.sim.matrix <- reshape2::melt(sim.matrix, id.vars = "organism")
  colnames(reshaped.sim.matrix) <- c("organism", "similarity","count")
  
  # df <- data.frame(sim = names(sim.matrix)[2:ncol(sim.matrix)],pvals = pvals.name)
  
  res <- ggplot2::ggplot(reshaped.sim.matrix, ggplot2::aes(x = similarity, y = count), order = FALSE) +
         ggplot2::geom_violin(scale = "width", ggplot2::aes(fill = similarity)) + 
         ggplot2::geom_jitter(width = 0.75, height = 0) + 
         ggplot2::theme_minimal() + ggplot2::labs(x = xlab, y = ylab, title = main) + 
         ggplot2::theme(legend.text = ggplot2::element_text(size = text.size)) + 
         ggplot2::theme(axis.title       = ggplot2::element_text(size = text.size,face = "bold"),
                        axis.text.y      = ggplot2::element_text(size = text.size,face = "bold"),
                        axis.text.x      = ggplot2::element_text(size = text.size,face = "bold"),
                        panel.background = ggplot2::element_blank(), 
                        plot.title       = ggplot2::element_text(size = text.size, colour = "black", face = "bold")) 
  # + ggplot2::geom_text(data = df, label = pvals)
  
  cat("Pairwise comparisons ( z-test ) :")
  cat("\n")
  pvals.res <- paste0("p = ",pvals)
  names(pvals.res) <- names(pvals)
  print(pvals.res)
  return (res)
}



# ifelse(!is.null(label.threshold),ggplot2::geom_text(ggplot2::aes(label = ifelse(count > label.threshold,as.character(organism),'')),hjust=0, vjust=0,size = 4), ggplot2::geom_text(ggplot2::aes(label = organism),hjust = 0, vjust = 0,size = 4))


