#' @title Plot Genome size vs. LTR transposon count
#' @description Genome size vs. LTR transposon count.
#' @param genome.matrix Genome matrix retruned by \code{\link{LTRpred.meta}}.
#' @param type type of LTR abundance. Options are: \code{type = "mass"} (total length of all TEs in Mbp), \code{type = "prop.mass"} (proportion of TEs within entire genome in %),
#' \code{type = "count"} (total number of TEs in genome), \code{type = "norm.count"} (total number of TEs in genome normalized by genome size in Mbp).
#' @param cor.method correlation analysis method. Either \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}.
#' @param smooth.method method to add a smoothed conditional mean (see \code{ggplot2::geom_smooth()}).
#' Users can choose from \code{auto}, \code{glm}, \code{lm}, or \code{loess}. If no smoothing function shall be applied,
#' users can specify \code{smooth.method = 'NULL'}. 
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
                                type                = "mass",
                                cor.method          = "spearman",
                                smooth.method       = "glm",
                                xlab                = "LTR transposon Abundance",
                                ylab                = "Genome size in Mega [bp]", 
                                main                = "Genome size vs. LTR transp. abundance",
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
    
    if (!is.element(type, c("mass", "prop.mass", "count", "norm.count")))
        stop ("Please specify: type = 'mass' for total length of all TEs in Mbp; type = 'prop.mass' for proportion of TEs within entire genome in %; type = 'count' for total number of TEs in genome; type = 'norm.count' for total number of TEs in genome normalized by genome size in Mbp.")
    
    if (cl.analysis) {
        if (is.null(sim.matrix))
            stop ("Please specify the 'sim.matrix' argument to be able to perform cluster analysis.")
        
        if (!is.null(smooth.method)) {
            if (!is.element(smooth.method, c("auto", "glm", "lm", "loess")))
                stop (
                    "Please specify a smooth.method that is supported by this function: 'auto', 'glm', 'lm', or 'loess'.",
                    call. = FALSE
                )
        }
        
        cl <-
            amap::Kmeans(
                x        = genome.matrix[, 2:ncol(genome.matrix)],
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
            cor(genome.matrix$nLTRs, genome.matrix$genome.size, method = cor.method)
       
         if (type == "norm.count")
            cor.value <-
            cor(genome.matrix$norm.nLTRs, genome.matrix$genome.size, method = cor.method)
        
        if (type == "mass")
            cor.value <-
                 cor(genome.matrix$totalMass, genome.matrix$genome.size, method = cor.method)
        
        if (type == "prop.mass")
            cor.value <-
            cor(genome.matrix$prop, genome.matrix$genome.size, method = cor.method)
        
        
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
                
                res <- res + ggplot2::geom_point() +
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
                ) + ggplot2::geom_text(
                    ggplot2::aes(label = organism),
                    hjust = 0,
                    vjust = 0,
                    size = label.size,
                    check_overlap = check.overlap
                )
            
            if (!is.null(smooth.method)) {
                res <- res + ggplot2::geom_smooth(method = smooth.method)
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
                res <- res + ggplot2::geom_point() +
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
                ) + ggrepel::geom_text_repel(ggplot2::aes(label = organism),
                                             size = 3,
                                             fontface = 'bold')
        }
        
        if (!is.null(smooth.method)) {
            res <- res + ggplot2::geom_smooth(method = smooth.method)
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
                
                res <- res + ggplot2::geom_point() +
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
                ) + ggplot2::geom_text(
                    ggplot2::aes(label = organism),
                    hjust = 0,
                    vjust = 0,
                    size = label.size,
                    check_overlap = check.overlap
                )
            
            if (!is.null(smooth.method)) {
                res <- res + ggplot2::geom_smooth(method = smooth.method)
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
                res <- res + ggplot2::geom_point() +
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
                ) + 
                ggrepel::geom_text_repel(ggplot2::aes(label = organism),size = 3,fontface = 'bold')
            
            if (!is.null(smooth.method)) {
                res <- res + ggplot2::geom_smooth(method = smooth.method)
            }
        }
        return (res)
    }
}

