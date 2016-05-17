#' @title Plot LTR Similarity vs. predicted LTR count
#' @description The LTR Similarity vs. predicted LTR count visualization
#' allows to study the variability of predicted LTRs in different genomes.
#' @param sim.matrix Similarity matrix retruned by \code{\link{LTRpred.meta}}.
#' @param genome.matrix Genome matrix retruned by \code{\link{LTRpred.meta}}.
#' @param type either absolute count value or normalized count value (normalized by genome size is shown): \code{type = "absolute"} or \code{type = "normalized"}. Default is \code{type = "absolute"}.
#' @param cl.analysis logical value indicating whether or not cluster analysis of the \code{sim.matrix} returned by \code{\link{LTRpred.meta}} shall be performed (see Details sections).
#' @param cl.centers number of expected clusters or a set of initial (distinct) cluster centres.
#' @param cl.method distance measure to perform cluster analysis. Options are \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"}, \code{"pearson"} , \code{"abspearson"} , \code{"abscorrelation"}, \code{"correlation"}, \code{"spearman"} or \code{"kendall"}.
#' @param cl.nstart If \code{cl.centers} is a number, number of random sets that shall be chosen.
#' @param cl.iter.max maximum number of iterations for cluster analysis.#' @param xlab x-axis label.
#' @param min.sim minimum similarity that was used with \code{\link{LTRpred.meta}}.
#' @param similarity.bin resolution of similarity binning that was used with \code{\link{LTRpred.meta}}.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param main main text.
#' @param text.size size of the labels in ggplot2 notation.
#' @author Hajk-Georg Drost
#' @details 
#' In case the \code{sim.matrix} includes more than 30 observations (genomes), then
#' a \code{z.test} is applied to compute pairwise differences between neighboring distributions.
#' 
#' In case the \code{sim.matrix} includes less than 30 observations (genomes), then
#' a \code{wilcox.test} is applied to compute pairwise differences between neighboring distributions.
#' 
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
                         genome.matrix       = NULL,
                         type                = "absolute",
                         cl.analysis         = FALSE,
                         cl.centers          = NULL,
                         cl.method           = "euclidean",
                         cl.nstart           = 100,
                         cl.iter.max         = 10000,
                         min.sim             = 70,
                         similarity.bin      = 2,
                         xlab                = "LTR % Similarity", 
                         ylab                = "Count", 
                         main                = "LTR % Similarity vs. Count",
                         text.size           = 18){
  
  
  if (!is.element(type, c("absolute", "normalized")))
    stop("Please choose a valid type: either 'absolute' or 'normalized'.", call. = FALSE)
  
  if ((type == "normalized") && (is.null(genome.matrix)))
    stop("Please specify the 'genome.matrix' argument when using type = 'normalized'.")
  
  similarity <- count <- NULL
  
  if (type == "normalized") {
    
    sim.matrix[ , 2:ncol(sim.matrix)] <- sim.matrix[ , 2:ncol(sim.matrix)] / genome.matrix$genome.size
    
  }
  
  names(sim.matrix) <- c("organism",levels(cut(100,rev(seq(100,min.sim,-similarity.bin)),
                                              include.lowest = TRUE,right = TRUE)))
  
  if (cl.analysis){
    
    cl <- amap::Kmeans(x        = sim.matrix[ , 2:ncol(sim.matrix)], 
                       centers  = cl.centers,
                       nstart   = cl.nstart, 
                       iter.max = cl.iter.max, 
                       method   = cl.method)
    
    colnames(sim.matrix)[1] <- "organism"
    reshaped.sim.matrix <- reshape2::melt(sim.matrix, id.vars = "organism")
    colnames(reshaped.sim.matrix) <- c("organism", "similarity","count")
    
    colors <- bcolor(max(cl$cluster))
    sim.matrix <- dplyr::mutate(sim.matrix,
                                cluster   = cl$cluster, 
                                cl.colors = colors[cl$cluster])
    
  } else {
    if (nrow(sim.matrix) > 30){
      cat("There are ",nrow(sim.matrix), " elements in the dataset.. \nA z.test is applied to compute pairwise distribution differences.")
      cat("\n")
      pvals <- pairwise.z.test(sim.matrix)
    }
    
    if (nrow(sim.matrix) < 30){
      cat("There are ", nrow(sim.matrix), " elements in the dataset.. \nA wilcox.test is applied to compute pairwise distribution differences.")
      cat("\n")
      pvals <- pairwise.wilcox.test(sim.matrix)
    }
    
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
    
    cat("Pairwise comparisons:")
    cat("\n")
    pvals.res <- paste0("p = ",round(pvals,digits = 2)," (",pvals.name,")")
    names(pvals.res) <- names(pvals)
    print(pvals.res)
    return (res)
  }
}



# ifelse(!is.null(label.threshold),ggplot2::geom_text(ggplot2::aes(label = ifelse(count > label.threshold,as.character(organism),'')),hjust=0, vjust=0,size = 4), ggplot2::geom_text(ggplot2::aes(label = organism),hjust = 0, vjust = 0,size = 4))


