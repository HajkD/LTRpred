#' @title Plot Genome size vs. LTR transposon count
#' @description Genome size vs. LTR transposon count.
#' @param genome.matrix Genome matrix retruned by \code{\link{LTRpred.meta}}.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param main main text.
#' @param text.size size of the labels in ggplot2 notation.
#' @author Hajk-Georg Drost
#' @examples 
#' \dontrun{
#' 
#' }
#' @export
PlotSizeCorrelation <- function(genome.matrix,
                                xlab      = "LTR transposon Count",
                                ylim      = "Genome size in bp", 
                                main      = "Genome size vs. LTR transp. count",
                                text.size = 18){
  
  colnames(genome.matrix)[1] <- "organism"
  res <- ggplot2::ggplot(genome.matrix, ggplot2::aes(x = nLTRs, y = genome.size)) + ggplot2::geom_point() + 
         ggplot2::theme_minimal() + 
         ggplot2::labs(x = xlab, y = ylim, title = main) + 
         ggplot2::theme(legend.text = ggplot2::element_text(size = text.size)) + 
         ggplot2::theme(axis.title  = ggplot2::element_text(size = text.size,face = "bold"),axis.text.y = ggplot2::element_text(size = text.size,face = "bold"),
                        axis.text.x = ggplot2::element_text(size = text.size,face = "bold"),
                        panel.background = ggplot2::element_blank(), 
                        plot.title = ggplot2::element_text(size = text.size, colour = "black", face = "bold"))
  
  return (res)
}

#          ggplot2::geom_text(ggplot2::aes(label = organism),hjust=0, vjust=0,size = 6) +


