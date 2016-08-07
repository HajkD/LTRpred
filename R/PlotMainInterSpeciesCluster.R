#' @title Plot inter species similarity between TEs (for the top 6 clusters)
#' @description After clustering TE sequences across different species with \code{\link{CLUSTpred}}
#' the similarity of TEs from different species for the top 6 clusters is being visualized.
#' @param cluster.file a cluster file in \code{*.uc} format (imported with e.g. \code{\link{filter.uc}} and \code{\link{read.uc}}). 
#' @param cluster numeric value denoting the cluster number of interest.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param legend.title legend title text.
#' @author Hajk-Georg Drost
#' @export

PlotMainInterSpeciesCluster <- function(cluster.file,
                                        xlab = "Species", 
                                        ylab = "% Similarity to Reference",
                                        legend.title = "Species"){
    # retrieve main clusters
    main.clusters <- as.numeric(names(head(sort(table(cluster.file$Clust_Cluster),decreasing = TRUE))))
    # generate PlotInterSpeciesCluster plots for the top 6 clusters
    p1 <- PlotInterSpeciesCluster(cluster.file, main.clusters[1], xlab = xlab, ylab = ylab, legend.title = legend.title)
    p2 <- PlotInterSpeciesCluster(cluster.file, main.clusters[2], xlab = xlab, ylab = ylab, legend.title = legend.title)
    p3 <- PlotInterSpeciesCluster(cluster.file, main.clusters[3], xlab = xlab, ylab = ylab, legend.title = legend.title)
    p4 <- PlotInterSpeciesCluster(cluster.file, main.clusters[4], xlab = xlab, ylab = ylab, legend.title = legend.title)
    p5 <- PlotInterSpeciesCluster(cluster.file, main.clusters[5], xlab = xlab, ylab = ylab, legend.title = legend.title)
    p6 <- PlotInterSpeciesCluster(cluster.file, main.clusters[6], xlab = xlab, ylab = ylab, legend.title = legend.title)
    
    cowplot::plot_grid(p1, p2, p3, p4, p5, p6)
}
