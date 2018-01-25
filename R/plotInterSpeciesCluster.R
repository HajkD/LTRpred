#' @title Plot inter species similarity between TEs (for a specific cluster)
#' @description After clustering TE sequences across different species with \code{\link[LTRpred]{CLUSTpred}}
#' the similarity of TEs from different species within one cluster can be visualized.
#' @param cluster.file a cluster file in \code{*.uc} format (imported with e.g. \code{\link[LTRpred]{filter.uc}} and \code{\link[LTRpred]{read.uc}}). 
#' @param cluster numeric value denoting the cluster number of interest.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param legend.title legend title text.
#' @author Hajk-Georg Drost
#' @export

plotInterSpeciesCluster <- function(cluster.file, 
                                    cluster,
                                    xlab = "Species", 
                                    ylab = "% Similarity to Reference",
                                    legend.title = "Species"){
    
    main.clusters <- as.numeric(names(sort(table(cluster.file$Clust_Cluster),decreasing = TRUE)))
    
    if (!is.element(cluster,main.clusters))
        stop("Cluster '",cluster,"' cannot be found in the cluster.file!", call. = FALSE)
    
    Clust_Cluster <- Species <- Clust_Perc_Ident <- n <- NULL
    
    # select all elements that belong to the same specific cluster
    cl <- dplyr::filter(cluster.file, Clust_Cluster == cluster)
    # define the species origin of the elements within this cluster
    cl <-
        dplyr::mutate(cl, Species = as.character(sapply(unlist(cl$Clust_Query), function(x)
            unlist(stringr::str_split(x, "_"))[1])))
    # count the total number elements for each species
    count <-
        as.data.frame(dplyr::summarise(dplyr::group_by(cl, Species), n = dplyr::n()))
    
    res <-
        ggplot2::ggplot(cl, ggplot2::aes(Species, Clust_Perc_Ident)) +
        ggplot2::geom_violin(ggplot2::aes(colour = Species), size = 2) +
        ggplot2::geom_point(ggplot2::aes(colour = Species)) +
        ggplot2::geom_jitter(ggplot2::aes(colour = Species),
                             position = ggplot2::position_jitter(0.3)) +
        ggplot2::labs(
            x = xlab,
            y = ylab,
            title = paste0("Cluster: ", cluster, " | Reference: ", unique(unlist(cl$Clust_Target)))
        ) +
        ggplot2::scale_fill_discrete(name = legend.title) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            title            = ggplot2::element_text(size = 10, face = "bold"),
            legend.title     = ggplot2::element_text(size = 18, face = "bold"),
            legend.text      = ggplot2::element_text(size = 18, face = "bold"),
            axis.title       = ggplot2::element_text(size = 18, face = "bold"),
            axis.text.y      = ggplot2::element_text(size = 18, face = "bold"),
            axis.text.x      = ggplot2::element_text(size = 18, face = "bold"),
            panel.background = ggplot2::element_blank(),
            strip.text.x     = ggplot2::element_text(
                size           = 18,
                colour         = "black",
                face           = "bold"
            )
        ) +
        ggplot2::ylim(c(min(cl$Clust_Perc_Ident) - 5, 100)) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 90,
            vjust = 1,
            hjust = 1
        )) +
        ggplot2::geom_text(
            data = count,
            ggplot2::aes(
                x = Species,
                y = min(cl$Clust_Perc_Ident) - 5,
                label = paste0("n = ", n),
                color = factor(Species)
            ),
            size = 6,
            fontface = "bold"
        )
    
    return(res)
}
