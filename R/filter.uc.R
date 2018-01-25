#' @title Filter for cluster members
#' @description Filter for cluster memebers in a \code{*.uc} format table.
#' @param cluster.file a \code{data.frame} \code{*.uc} (USEARCH cluster) format (imported with \code{\link{read.uc}}).
#' @author Hajk-Georg Drost
#' @export

filter.uc <- function(cluster.file){
     
    Type <- Cluster <- Query <- Target <- Perc_Ident <- NULL
     if (nrow(cluster.file) > 0) {
         cluster.file.type.H <-
             dplyr::filter(cluster.file, Type == "H")
         cluster.file.type.H <-
             dplyr::select(cluster.file.type.H,
                           Cluster,
                           Query,
                           Target,
                           Perc_Ident)
         cluster.file.type.H <-
             dplyr::mutate(cluster.file.type.H, Perc_Ident = as.numeric(Perc_Ident))
         names(cluster.file.type.H) <-
             paste0("Clust_", names(cluster.file.type.H))
     } else {
             stop("Your input cluster is empty.", call. = FALSE)
     }
     
     return(cluster.file.type.H)
 }
