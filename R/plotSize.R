#' @title Genome size vs. LTR transposon content
#' @description Compute and visualize the correlation between genome size vs. LTR transposon content.
#' @param genome.matrix Genome matrix retruned by \code{\link[LTRpred]{LTRpred.meta}}.
#' @param type type of LTR abundance. Options are: \code{type = "mass"} (total length of all TEs in Mbp), \code{type = "prop.mass"} (proportion of TEs within entire genome in %),
#' \code{type = "count"} (total number of TEs in genome), \code{type = "norm.count"} (total number of TEs in genome normalized by genome size in Mbp).
#' @param cor.method correlation analysis method. Either \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}.
#' @param label.organism shall organism names be drawn in the figure?
#' @param smooth.method method to add a smoothed conditional mean (see \code{ggplot2::geom_smooth()}).
#' Users can choose from \code{auto}, \code{glm}, \code{lm}, or \code{loess}. If no smoothing function shall be applied,
#' users can specify \code{smooth.method = 'NULL'}.
#' @param se shall standard error of the smoothing model (grey area) be drawn or not. 
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param main main text.
#' @param alpha transparency of confidence interval of the smoothing function: between \[ 0,1 \].
#' @param arrow_lab whether or not data points shall be labelled with arrows (ggrepel). 
#' @param cl.analysis logical value indicating whether or not cluster analysis of the \code{sim.matrix} returned by \code{\link[LTRpred]{LTRpred.meta}} shall be performed (see Details sections).
#' @param cl.centers number of expected clusters or a set of initial (distinct) cluster centres.
#' @param cl.method distance measure to perform cluster analysis. Options are \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"}, \code{"pearson"} , \code{"abspearson"} , \code{"abscorrelation"}, \code{"correlation"}, \code{"spearman"} or \code{"kendall"}.
#' @param cl.nstart If \code{cl.centers} is a number, number of random sets that shall be chosen.
#' @param cl.iter.max maximum number of iterations for cluster analysis.
#' @param sim.matrix Similarity matrix retruned by \code{\link[LTRpred]{LTRpred.meta}}.
#' @param text.size size of the labels in ggplot2 notation.
#' @param label.size size of the dot-labels.
#' @param check.overlap shall overlaps of labels be avoided or not.
#' @author Hajk-Georg Drost
#' @details 
#' In the first step, users should have generated LTR retrotransposon annotations for several species using the
#' \code{\link[LTRpred]{LTRpred.meta}} function. The output of \code{\link[LTRpred]{LTRpred.meta}} is a folder storing the annotation files
#' for all species of input species. In addition, \code{\link[LTRpred]{LTRpred.meta}} generates two files: \code{*_SimilarityMatrix} and \code{*_GenomeInfo} .
#' @examples 
#' GenomeMatrix <- read.csv(system.file("GenomeMatrix.csv",package = "LTRpred"), sep = ";")
#' GenomeMatrix[ , 1] <- sapply(GenomeMatrix[ , 1], 
#'                              function(x) unlist(stringr::str_split(x,"_"))[1])
#'                              
#' plotSize(GenomeMatrix)
#' @export

plotSize <- function(genome.matrix,
                                type                = "mass",
                                cor.method          = "spearman",
                                label.organism      = TRUE,
                                smooth.method       = "glm",
                                se                  = TRUE,
                                xlab                = "LTR retrotransposon content in Mega [bp]",
                                ylab                = "Genome size in Mega [bp]", 
                                main                = "",
                                alpha               = 0.3,
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
  
  
    nLTRs <- genome.size <- totalMass <- cl.colors <- prop <- organism <- norm.nLTRs <- NULL
    
    if (!is.element(type, c("mass", "prop.mass", "count", "norm.count")))
        stop("Please specify: type = 'mass' for total length of all TEs in Mbp; type = 'prop.mass' for proportion of TEs within entire genome in %; type = 'count' for total number of TEs in genome; type = 'norm.count' for total number of TEs in genome normalized by genome size in Mbp.")
    
    if (cl.analysis) {
        if (is.null(sim.matrix))
            stop("Please specify the 'sim.matrix' argument to be able to perform cluster analysis.")
        
        if (!is.null(smooth.method)) {
            if (!is.element(smooth.method, c("auto", "glm", "lm", "loess")))
                stop(
                    "Please specify a smooth.method that is supported by this function: 'auto', 'glm', 'lm', or 'loess'.",
                    call. = FALSE
                )
        }
        
        cl <-
            amap::Kmeans(
                x        = LTRpred::SimMatAbundance(sim.matrix),
                centers  = cl.centers,
                nstart   = cl.nstart,
                iter.max = cl.iter.max,
                method   = cl.method
            )
        
        colors <- bcolor(max(cl$cluster))
        genome.matrix <-
            dplyr::mutate(genome.matrix, cl.colors = colors[cl$cluster])
        
        if (type == "count")
        cor.value <-
            stats::cor(genome.matrix$nLTRs, genome.matrix$genome.size, method = cor.method)
       
         if (type == "norm.count")
            cor.value <-
            stats::cor(genome.matrix$norm.nLTRs, genome.matrix$genome.size, method = cor.method)
        
        if (type == "mass")
            cor.value <-
            stats::cor(genome.matrix$totalMass, genome.matrix$genome.size, method = cor.method)
        
        if (type == "prop.mass")
            cor.value <-
            stats::cor(genome.matrix$prop, genome.matrix$genome.size, method = cor.method)
        
        
        colnames(genome.matrix)[1] <- "organism"
        
        if (!arrow_lab) {
            
            if (type == "count")
            res <-
                ggplot2::ggplot(genome.matrix,
                                ggplot2::aes(
                                    x = nLTRs,
                                    y = genome.size,
                                    colour = cl.colors
                                ))  
            
            if (type == "norm.count")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(
                                        x = norm.nLTRs,
                                        y = genome.size,
                                        colour = cl.colors
                                    ))  
            
            if (type == "mass")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(
                                        x = totalMass,
                                        y = genome.size,
                                        colour = cl.colors
                                    ))    
            if (type == "prop.mass")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(
                                        x = prop * 100,
                                        y = genome.size,
                                        colour = cl.colors
                                    )) 
                
                res <- res + ggplot2::geom_point(size = 5) +
                ggplot2::theme_minimal() +
                ggplot2::labs(
                    x = xlab,
                    y = ylab,
                    title = paste0(
                        main,
                        " [ Corr ( ",
                        cor.method,
                        " ) = ",
                        round(cor.value, digits = 2),
                        " ]"
                    )
                ) +
                ggplot2::theme(legend.text = ggplot2::element_text(size = text.size)) +
                ggplot2::theme(
                    axis.title  = ggplot2::element_text(size = text.size, face = "bold"),
                    axis.text.y = ggplot2::element_text(size = text.size, face = "bold"),
                    axis.text.x = ggplot2::element_text(size = text.size, face = "bold"),
                    panel.background = ggplot2::element_blank(),
                    plot.title = ggplot2::element_text(
                        size = text.size,
                        colour = "black",
                        face = "bold"
                    )
                ) 
                
                if (label.organism) {
                    res <- res + ggplot2::geom_text(
                        ggplot2::aes(label = organism),
                        hjust = 0,
                        vjust = 0,
                        size = label.size,
                        check_overlap = check.overlap
                    ) 
                }
                res <- res +
                    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
                    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6))
            
            if (!is.null(smooth.method)) {
                res <- res + ggplot2::geom_smooth(method = smooth.method, se = se, alpha = alpha)
            }
        }
        
        if (arrow_lab) {
            
            if (type == "count")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(x = nLTRs,
                                                 y = genome.size,
                                                 colour = cl.colors)) 
            
            if (type == "norm.count")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(x = norm.nLTRs,
                                                 y = genome.size,
                                                 colour = cl.colors)) 
            
            if (type == "mass")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(
                                        x = totalMass,
                                        y = genome.size,
                                        colour = cl.colors
                                    ))    
            if (type == "prop.mass")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(
                                        x = prop * 100,
                                        y = genome.size,
                                        colour = cl.colors
                                    )) 
                res <- res + ggplot2::geom_point(size = 5) +
                ggplot2::theme_minimal() +
                ggplot2::labs(
                    x = xlab,
                    y = ylab,
                    title = paste0(
                        main,
                        " [ Corr ( ",
                        cor.method,
                        " ) = ",
                        round(cor.value, digits = 2),
                        " ]"
                    )
                ) +
                ggplot2::theme(legend.text = ggplot2::element_text(size = text.size)) +
                ggplot2::theme(
                    axis.title  = ggplot2::element_text(size = text.size, face = "bold"),
                    axis.text.y = ggplot2::element_text(size = text.size, face = "bold"),
                    axis.text.x = ggplot2::element_text(size = text.size, face = "bold"),
                    panel.background = ggplot2::element_blank(),
                    plot.title = ggplot2::element_text(
                        size = text.size,
                        colour = "black",
                        face = "bold"
                    )
                ) 
                
                if (label.organism) {
                    res <- res + ggrepel::geom_text_repel(ggplot2::aes(label = organism),
                                                          size = label.size,
                                                          fontface = 'bold')
                }
                res <- res +
                    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
                    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6))
        }
        
        if (!is.null(smooth.method)) {
            res <- res + ggplot2::geom_smooth(method = smooth.method, se = se, alpha = alpha)
        }
        
        print(res)
        return (cl)
    } else {
        
        if (type == "count")
            cor.value <-
                stats::cor(genome.matrix$nLTRs, genome.matrix$genome.size, method = cor.method)
        if (type == "norm.count")
            cor.value <-
                stats::cor(genome.matrix$norm.nLTRs, genome.matrix$genome.size, method = cor.method)
        if (type == "mass")
            cor.value <-
                stats::cor(genome.matrix$totalMass, genome.matrix$genome.size, method = cor.method)
        if (type == "prop.mass")
            cor.value <-
                stats::cor(genome.matrix$prop, genome.matrix$genome.size, method = cor.method)
        
        colnames(genome.matrix)[1] <- "organism"
        
        if (!arrow_lab) {
            
            if (type == "count")
               res <-
                   ggplot2::ggplot(genome.matrix,
                                ggplot2::aes(
                                    x = nLTRs,
                                    y = genome.size
                                ))  
            
            if (type == "norm.count")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(
                                        x = norm.nLTRs,
                                        y = genome.size
                                    ))  
            
            if (type == "mass")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(
                                        x = totalMass,
                                        y = genome.size
                                    ))    
            if (type == "prop.mass")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(
                                        x = prop * 100,
                                        y = genome.size
                                    )) 
                
                res <- res + ggplot2::geom_point(size = 5) +
                ggplot2::theme_minimal() +
                ggplot2::labs(
                    x = xlab,
                    y = ylab,
                    title = paste0(
                        main,
                        " [ Corr ( ",
                        cor.method,
                        " ) = ",
                        round(cor.value, digits = 2),
                        " ]"
                    )
                ) +
                ggplot2::theme(legend.text = ggplot2::element_text(size = text.size)) +
                ggplot2::theme(
                    axis.title  = ggplot2::element_text(size = text.size, face = "bold"),
                    axis.text.y = ggplot2::element_text(size = text.size, face = "bold"),
                    axis.text.x = ggplot2::element_text(size = text.size, face = "bold"),
                    panel.background = ggplot2::element_blank(),
                    plot.title = ggplot2::element_text(
                        size = text.size,
                        colour = "black",
                        face = "bold"
                    )
                ) 
                
                if (label.organism) {
                    res <- res + ggplot2::geom_text(
                        ggplot2::aes(label = organism),
                        hjust = 0,
                        vjust = 0,
                        size = label.size,
                        check_overlap = check.overlap
                    )
                }
                
                res <- res +
                    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
                    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6))
            
            if (!is.null(smooth.method)) {
                res <- res + ggplot2::geom_smooth(method = smooth.method, se = se, alpha = alpha)
            }
        }
        
        if (arrow_lab) {
            if (type == "count")
               res <-
                ggplot2::ggplot(genome.matrix,
                                ggplot2::aes(
                                    x = nLTRs,
                                    y = genome.size
                                ))  
            if (type == "norm.count")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(
                                        x = norm.nLTRs,
                                        y = genome.size
                                    )) 
            if (type == "mass")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(
                                        x = totalMass,
                                        y = genome.size
                                    ))    
            if (type == "prop.mass")
                res <-
                    ggplot2::ggplot(genome.matrix,
                                    ggplot2::aes(
                                        x = prop * 100,
                                        y = genome.size
                                    ))    
                res <- res + ggplot2::geom_point(size = 5) +
                ggplot2::theme_minimal() +
                ggplot2::labs(
                    x = xlab,
                    y = ylab,
                    title = paste0(
                        main,
                        " [ Corr ( ",
                        cor.method,
                        " ) = ",
                        round(cor.value, digits = 2),
                        " ]"
                    )
                ) +
                ggplot2::theme(legend.text = ggplot2::element_text(size = text.size)) +
                ggplot2::theme(
                    axis.title  = ggplot2::element_text(size = text.size, face = "bold"),
                    axis.text.y = ggplot2::element_text(size = text.size, face = "bold"),
                    axis.text.x = ggplot2::element_text(size = text.size, face = "bold"),
                    panel.background = ggplot2::element_blank(),
                    plot.title = ggplot2::element_text(
                        size = text.size,
                        colour = "black",
                        face = "bold"
                    )
                ) 
                
                if (label.organism) {
                    res <- res + ggrepel::geom_text_repel(ggplot2::aes(label = organism), size = label.size,fontface = 'bold')
                }
                
                    res <- res + ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
                    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6))
            
            if (!is.null(smooth.method)) {
                res <- res + ggplot2::geom_smooth(method = smooth.method, se = se, alpha = alpha)
            }
        }
        return (res)
    }
}

