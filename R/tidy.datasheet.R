#' @title Select most important columns of \code{LTRpred} output for further analytics
#' @description This function simply selects the most important columns of \code{LTRpred} output for further analytics.
#' @param LTRpred.tbl a \code{data.frame} storing the result (DataSheet) of \code{\link{LTRpred}}.
#' @author Hajk-Georg Drost
#' @export

tidy.datasheet <- function(LTRpred.tbl){
    
    ID <- chromosome <- start <- end <- strand <- dfam_target_name <- dfam_acc <- width <- NULL
    ltr_similarity <- orfs <- protein_domain <- Clust_Cluster <- Clust_Target <- Clust_Perc_Ident <- NULL
    Clust_cn <- cn_3ltr <- cn_5ltr <- lLTR_start <- `PBS/tRNA_edist` <- PPT_length <- PBS_length <- dfam_target_description <- NULL
    orf.id <- TE_CG_abs <- N_5ltr_abs <- NULL
    
    
    return(
        dplyr::select(
            LTRpred.tbl,
            ID,
            chromosome,
            start,
            end,
            strand,
            dfam_target_name,
            dfam_acc,
            width,
            ltr_similarity,
            orfs,
            protein_domain,
            Clust_Cluster,
            Clust_Target,
            Clust_Perc_Ident,
            Clust_cn,
            cn_3ltr,
            cn_5ltr,
            lLTR_start:`PBS/tRNA_edist`,
            PPT_length,
            PBS_length,
            dfam_target_description,
            orf.id,
            TE_CG_abs:N_5ltr_abs
        )
    )
    
    
}
