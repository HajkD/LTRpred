#' @title Plot correlation between Cluster Number and any other variable
#' @description Visualize the correlation between the Cluster Number generated 
#' by \code{\link{LTRpred}} and any other variable 
#' and color each putative LTR transposon by it's e.g. superfamily name.
#' @param LTRpred.tbl \code{data.frame} returned by \code{\link{LTRpred}}.
#' @param y column name in \code{LTRpred.tbl} (lazy evaluation) that shall be correlated with \code{cn_3ltr} (copy number of predicted 3' LTR).
#' @param color.by color points by a column name in \code{LTRpred.tbl} (lazy evaluation).
#' @param quality.filter shall a quality filter be applied to reduce false positives?
#' @param sim filter minimum sequence similarity between the LTRs.
#' @param n.orfs filter for minimum number of orfs.
#' @param chop.names shall superfamily names be chopped?
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param main main text.
#' @param legend.title legend text.
#' @author Hajk-Georg Drost
#' @export

PlotCluster <-
    function(LTRpred.tbl,
             y = "CG_3ltr_rel",
             color.by = "fam",
             quality.filter = FALSE,
             sim = 70,
             n.orfs = 1,
             chop.names = TRUE,
             xlab = "Cluster Number",
             ylab = "",
             main = "",
             legend.title = "Protein Domain") {
        
        if (!is.element(color.by, c("fam", "dfam_target_name")))
            stop("Please specify either 'color.by = 'dfam_target_name' or 'color.by = 'fam' (in case 'chop.names = TRUE').")
        
        if (chop.names && (color.by != "fam"))
            stop("When specifying 'chop.names = FALSE' please also specify 'color.by = 'dfam_target_name'.'")
        
        ltr_similarity <- PBS_start <- protein_domain <- orfs <- dfam_target_name <- Clust_Cluster <- NULL 
        
        # chop superfamily names
        fam_replace <- function(x){
            
            if (stringr::str_detect(x,"Copia")) {
                return(stringr::str_replace(x,x,"Copia"))
            }
            if (stringr::str_detect(x,"Gypsy")) {
                return(stringr::str_replace(x,x,"Gypsy"))
            }
            if (stringr::str_detect(x,"DNA")) {
                return(stringr::str_replace(x,x,"DNA"))
            }
            
            return(x)
        }
        
        if (quality.filter) {
            LTRpred.tbl <- dplyr::filter(LTRpred.tbl,
                                         ltr_similarity >= sim,
                                         (!is.na(PBS_start)) |
                                             (!is.na(protein_domain)),
                                         orfs >= n.orfs)
        }
        
        if (chop.names)
            LTRpred.tbl <- dplyr::mutate(LTRpred.tbl, fam = unlist(sapply(dfam_target_name,fam_replace)))
        
        
        cl_names <- names(table(LTRpred.tbl$Clust_Cluster))
        cl_names <- cl_names[which(cl_names != "unique")]

        cluster.df <- do.call(rbind,lapply(cl_names, function(x) {
            
            cluster.members(LTRpred.tbl, x)
            
        }))
        
        cluster.df.final <- rbind(cluster.df, dplyr::filter(LTRpred.tbl, Clust_Cluster == "unique"))
        
        p <- ggplot2::ggplot(cluster.df.final, ggplot2::aes_string(x = "Clust_Cluster", y = y, color = color.by)) +
            ggplot2::geom_point() +
            ggplot2::theme_minimal() + 
            ggplot2::labs(x = xlab, y = ylab, title = main) +
            ggplot2::scale_fill_discrete(name = legend.title) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                title            = ggplot2::element_text(size = 18, face = "bold"),
                legend.title     = ggplot2::element_text(size = 18, face = "bold"),
                legend.text      = ggplot2::element_text(size = 18, face = "bold"),
                axis.title       = ggplot2::element_text(size = 18, face = "bold"),
                axis.text.y      = ggplot2::element_text(size = 18, face = "bold"),
                axis.text.x      = ggplot2::element_text(size = 14, face = "bold"),
                panel.background = ggplot2::element_blank(),
                strip.text.x     = ggplot2::element_text(
                    size           = 18,
                    colour         = "black",
                    face           = "bold"
                )
            ) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
        
        
        return(p)
    }














