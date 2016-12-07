#' @title Plot inter species similarity between TEs (for the top 6 clusters)
#' @description After clustering TE sequences across different species with \code{\link{CLUSTpred}}
#' the similarity of TEs from different species for the top 6 clusters is being visualized.
#' @param cluster.file a cluster file in \code{*.uc} format (imported with e.g. \code{\link{filter.uc}} and \code{\link{read.uc}}).
#' @param top.n number of top clustern that shall be visualized.
#' @param decreasing shall clusters be shown by size in decreasing order? 
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param legend.title legend title text.
#' @author Hajk-Georg Drost
#' @export

PlotMainInterSpeciesCluster <- function(cluster.file,
                                        top.n = 6,
                                        decreasing = TRUE,
                                        xlab = "Species", 
                                        ylab = "% Similarity to Reference",
                                        legend.title = "Species"){
    # retrieve main clusters
    main.clusters <-
        as.numeric(names(sort(
            table(cluster.file$Clust_Cluster), decreasing = decreasing
        )))
    
    if (length(main.clusters) < top.n)
        stop(
            "Number of clusters [",
            length(main.clusters),
            "] is smaller than number of selected elements 'top.n = ",
            top.n,
            "'.",
            call. = FALSE
        )
    
    nf <- grDevices::n2mfrow(top.n)
    pl <-
        lapply(seq_len(top.n), function(x)
            PlotInterSpeciesCluster(
                cluster.file,
                main.clusters[x],
                xlab = xlab,
                ylab = ylab,
                legend.title = legend.title
            ))
    gridExtra::marrangeGrob(pl, nrow = nf[1], ncol = nf[2])
}
