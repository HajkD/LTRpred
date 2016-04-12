#' @title Plot Genome size vs. LTR transposon count
#' @description Genome size vs. LTR transposon count.
#' @param genome.matrix Genome matrix retruned by \code{\link{LTRpred.meta}}.
#' @param cor.method correlation analysis method. Either \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}. 
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param main main text.
#' @param arrow_lab whether or not data points shall be labelled with arrows (ggrepel). 
#' @param cl.analysis logical value indicating whether or not cluster analysis of the \code{sim.matrix} returned by \code{\link{LTRpred.meta}} shall be performed (see Details sections).
#' @param cl.centers number of expected clusters or a set of initial (distinct) cluster centres.
#' @param cl.method distance measure to perform cluster analysis. Options are \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"}, \code{"pearson"} , \code{"abspearson"} , \code{"abscorrelation"}, \code{"correlation"}, \code{"spearman"} or \code{"kendall"}.
#' @param cl.nstart If \code{cl.centers} is a number, number of random sets that shall be chosen.
#' @param cl.iter.max maximum number of iterations for cluster analysis.
#' @param sim.matrix Similarity matrix retruned by \code{\link{LTRpred.meta}}.
#' @param text.size size of the labels in ggplot2 notation.
#' @param label.size size of the dot-labels.
#' @param check.overlap shall overlaps of labels be avoided or not.
#' @author Hajk-Georg Drost
#' @examples 
#' GenomeMatrix <- read.csv(system.file("GenomeMatrix.csv",package = "LTRpred"), sep = ";")
#' GenomeMatrix[ , 1] <- sapply(GenomeMatrix[ , 1], 
#'                              function(x) unlist(stringr::str_split(x,"_"))[1])
#'                              
#' PlotSizeCorrelation(GenomeMatrix)
#' @export

PlotSizeCorrelation <- function(genome.matrix,
                                cor.method          = "spearman",
                                xlab                = "LTR transposon Count",
                                ylab                = "Genome size in bp", 
                                main                = "Genome size vs. LTR transp. count",
                                arrow_lab           = FALSE,
                                cl.analysis         = FALSE,
                                cl.centers          = NULL,
                                cl.method           = "euclidean",
                                cl.nstart           = 100,
                                cl.iter.max         = 10000,
                                sim.matrix          = NULL,
                                text.size           = 18,
                                label.size          = 3,
                                check.overlap       = TRUE){
  
  
  nLTRs <- genome.size <- cl.colors <- organism <- NULL
  
  if (cl.analysis){
    
    if (is.null(sim.matrix))
      stop ("Please specify the 'sim.matrix' argument to be able to perform cluster analysis.")
    
    cl <- amap::Kmeans(x        = sim.matrix[ , 2:ncol(sim.matrix)], 
                       centers  = cl.centers,
                       nstart   = cl.nstart, 
                       iter.max = cl.iter.max, 
                       method   = cl.method)
    
    colors <- bcolor(max(cl$cluster))
    genome.matrix <- dplyr::mutate(genome.matrix, cl.colors = colors[cl$cluster])
    
    cor.value <- cor(genome.matrix$nLTRs,genome.matrix$genome.size, method = cor.method)
    
    colnames(genome.matrix)[1] <- "organism"
    
    if (!arrow_lab){
      res <- ggplot2::ggplot(genome.matrix, ggplot2::aes(x = nLTRs, y = genome.size, colour = cl.colors)) + ggplot2::geom_point() + 
        ggplot2::theme_minimal() + 
        ggplot2::labs(x = xlab, y = ylab, title = paste0(main, " [ Corr ( ",cor.method," ) = ",round(cor.value, digits = 2), " ]")) + 
        ggplot2::theme(legend.text = ggplot2::element_text(size = text.size)) + 
        ggplot2::theme(axis.title  = ggplot2::element_text(size = text.size,face = "bold"),axis.text.y = ggplot2::element_text(size = text.size,face = "bold"),
                       axis.text.x = ggplot2::element_text(size = text.size,face = "bold"),
                       panel.background = ggplot2::element_blank(), 
                       plot.title = ggplot2::element_text(size = text.size, colour = "black", face = "bold")) + ggplot2::geom_text(ggplot2::aes(label = organism), hjust = 0, vjust = 0, 
                                                                                                                                   size = label.size, check_overlap = check.overlap)
    }
    
    if (arrow_lab){
      
      res <- ggplot2::ggplot(genome.matrix, ggplot2::aes(x = nLTRs, y = genome.size, colour = cl.colors)) + ggplot2::geom_point() + 
        ggplot2::theme_minimal() + 
        ggplot2::labs(x = xlab, y = ylab, title = paste0(main, " [ Corr ( ",cor.method," ) = ",round(cor.value, digits = 2), " ]")) + 
        ggplot2::theme(legend.text = ggplot2::element_text(size = text.size)) + 
        ggplot2::theme(axis.title  = ggplot2::element_text(size = text.size,face = "bold"),axis.text.y = ggplot2::element_text(size = text.size,face = "bold"),
                       axis.text.x = ggplot2::element_text(size = text.size,face = "bold"),
                       panel.background = ggplot2::element_blank(), 
                       plot.title = ggplot2::element_text(size = text.size, colour = "black", face = "bold")) + ggrepel::geom_text_repel(ggplot2::aes(label = organism), size = 2,fontface = 'bold')
    }
    
    print(res)
    return (cl)
  } else {
      
    cor.value <- stats::cor(genome.matrix$nLTRs, genome.matrix$genome.size, method = cor.method)
    
    colnames(genome.matrix)[1] <- "organism"
    
    if (!arrow_lab){
      
      res <- ggplot2::ggplot(genome.matrix, ggplot2::aes(x = nLTRs, y = genome.size)) + ggplot2::geom_point() + 
        ggplot2::theme_minimal() + 
        ggplot2::labs(x = xlab, y = ylab, title = paste0(main, " [ Corr ( ",cor.method," ) = ",round(cor.value, digits = 2), " ]")) + 
        ggplot2::theme(legend.text = ggplot2::element_text(size = text.size)) + 
        ggplot2::theme(axis.title  = ggplot2::element_text(size = text.size,face = "bold"),axis.text.y = ggplot2::element_text(size = text.size,face = "bold"),
                       axis.text.x = ggplot2::element_text(size = text.size,face = "bold"),
                       panel.background = ggplot2::element_blank(), 
                       plot.title = ggplot2::element_text(size = text.size, colour = "black", face = "bold")) + ggplot2::geom_text(ggplot2::aes(label = organism), hjust = 0, vjust = 0, 
                                                                                                                                   size = label.size, check_overlap = check.overlap)
    }
    
    if (arrow_lab){
      res <- ggplot2::ggplot(genome.matrix, ggplot2::aes(x = nLTRs, y = genome.size)) + ggplot2::geom_point() + 
        ggplot2::theme_minimal() + 
        ggplot2::labs(x = xlab, y = ylab, title = paste0(main, " [ Corr ( ",cor.method," ) = ",round(cor.value, digits = 2), " ]")) + 
        ggplot2::theme(legend.text = ggplot2::element_text(size = text.size)) + ggrepel::geom_text_repel(ggplot2::aes(label = organism), size = 2, fontface = 'bold')
    }
    
    return (res)
    
    }
}

