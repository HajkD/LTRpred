#' @title Select most important columns of \code{LTRpred} output for further analytics
#' @description This function simply selects the most important columns of \code{LTRpred} output for further analytics.
#' @param LTRpred.tbl a \code{data.frame} storing the result (DataSheet) of \code{\link{LTRpred}}.
#' @author Hajk-Georg Drost
#' @export

tidy.datasheet <- function(LTRpred.tbl){
    
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
            cn_3ltr,
            cn_5ltr,
            lLTR_start:`PBS/tRNA_edist`,
            PPT_length,
            PBS_length,
            dfam_target_description,
            orf.id,
            TE_CG_abs:CHH_5ltr_rel
        )
    )
    
    
}
