#' @title Import \code{LTRpred} DataSheet
#' @description The \code{*_LTRpred_DataSheet.tsv} file generated by \code{\link{LTRpred}} stores the features of all predicted LTR transposons in a table. This function
#' imports this \code{\link{LTRpred}} output table.
#' @param data.sheet path to the \code{*_LTRpred_DataSheet.tsv} file.
#' @author Hajk-Georg Drost
#' @examples
#' # example prediction file generated by LTRpred 
#' pred.file <- system.file("Hsapiens_ChrY_LTRpred_DataSheet.tsv", package = "LTRpred")
#' # read LTRpred generated prediction file (data sheet)
#' pred <- read.ltrpred(pred.file)
#' 
#' # or arrange by ltr_similarity
#' dplyr::arrange(tidy.datasheet(pred), dplyr::desc(ltr_similarity))
#' @export
   
read.ltrpred <- function(data.sheet){
    
    if(!file.exists(data.sheet))
        stop("The file ", data.sheet, " does not seem to exist. Please provide a valid path to a file named *_LTRpred_DataSheet.tsv.", call. = FALSE)
    
    if (file.info(data.sheet)$size == 0 ||
        is.na(file.info(data.sheet)$size == 0)) {
        cat("File ",
            data.sheet,
            " is empty and therefore is not being processed.")
        cat("\n")
        cat("\n")
        return(NULL)
    }
    
    get_header <- readr::read_lines(data.sheet, n_max = 1)
    get_header_clean <- unlist(stringr::str_split(get_header, "\t"))
    
    if (all(c("ltr_age_mya", "ltr_evo_distance") %in% get_header_clean)) {
        pred <- readr::read_delim(data.sheet, col_types = readr::cols(
            "ID" = readr::col_character(),
            "dfam_target_name" = readr::col_character(),
            "ltr_similarity" = readr::col_double(),
            "ltr_age_mya" = readr::col_double(),
            "ltr_evo_distance" = readr::col_double(),
            "similarity" = readr::col_character(),
            "protein_domain" = readr::col_character(),
            "orfs" = readr::col_integer(),
            "chromosome" = readr::col_character(),
            "start" = readr::col_integer(),
            "end" = readr::col_integer(),
            "strand" = readr::col_character(),
            "width" = readr::col_integer(),
            "annotation" = readr::col_character(),
            "pred_tool" = readr::col_character(),
            "frame" = readr::col_character(),
            "score" = readr::col_character(),
            "lLTR_start" = readr::col_integer(),
            "lLTR_end" = readr::col_integer(),
            "lLTR_length" = readr::col_integer(),
            "rLTR_start" = readr::col_integer(),
            "rLTR_end" = readr::col_integer(),
            "rLTR_length" = readr::col_integer(),
            "lTSD_start" = readr::col_integer(),
            "lTSD_end" = readr::col_integer(),
            "lTSD_motif" = readr::col_character(),
            "rTSD_start" = readr::col_integer(),
            "rTSD_end" = readr::col_integer(),
            "rTSD_motif" = readr::col_character(),
            "PPT_start" = readr::col_integer(),
            "PPT_end" = readr::col_integer(),
            "PPT_motif" = readr::col_character(),
            "PPT_strand" = readr::col_character(),
            "PPT_offset" = readr::col_integer(),
            "PBS_start" = readr::col_integer(),
            "PBS_end" = readr::col_integer(),
            "PBS_strand" = readr::col_character(),
            "tRNA" = readr::col_character(),
            "tRNA_motif" = readr::col_character(),
            "PBS_offset" = readr::col_integer(),
            "tRNA_offset" = readr::col_integer(),
            "PBS/tRNA_edist" = readr::col_integer(),
            "orf.id" = readr::col_character(),
            "repeat_region_length" = readr::col_integer(),
            "PPT_length" = readr::col_integer(),
            "PBS_length" = readr::col_integer(),
            "dfam_acc" = readr::col_character(),
            "dfam_bits" = readr::col_double(),
            "dfam_e_value" = readr::col_double(),
            "dfam_bias" = readr::col_double(),
            "dfam_hmm-st" = readr::col_double(),
            "dfam_hmm-en" = readr::col_double(),
            "dfam_strand" = readr::col_character(),
            "dfam_ali-st" = readr::col_double(),
            "dfam_ali-en" = readr::col_double(),
            "dfam_env-st" = readr::col_double(),
            "dfam_env-en" = readr::col_double(),
            "dfam_modlen" = readr::col_double(),
            "dfam_target_description" = readr::col_character(),
            "Clust_Cluster" = readr::col_character(),
            "Clust_Target" = readr::col_character(),
            "Clust_Perc_Ident" = readr::col_double(),
            "Clust_cn" = readr::col_integer(),
            "TE_CG_abs" = readr::col_double(),
            "TE_CG_rel" = readr::col_double(),
            "TE_CHG_abs" = readr::col_double(),
            "TE_CHG_rel" = readr::col_double(),
            "TE_CHH_abs" = readr::col_double(),
            "TE_CHH_rel" = readr::col_double(),
            "TE_CCG_abs" = readr::col_double(),
            "TE_CCG_rel" = readr::col_double(),
            "TE_N_abs" = readr::col_double(),
            "CG_3ltr_abs" = readr::col_double(),
            "CG_3ltr_rel" = readr::col_double(),
            "CHG_3ltr_abs" = readr::col_double(),
            "CHG_3ltr_rel" = readr::col_double(),
            "CHH_3ltr_abs" = readr::col_double(),
            "CHH_3ltr_rel" = readr::col_double(),
            "CCG_3ltr_abs" = readr::col_double(),
            "CCG_3ltr_rel" = readr::col_double(),
            "N_3ltr_abs" = readr::col_double(),
            "CG_5ltr_abs" = readr::col_double(),
            "CG_5ltr_rel" = readr::col_double(),
            "CHG_5ltr_abs" = readr::col_double(),
            "CHG_5ltr_rel" = readr::col_double(),
            "CHH_5ltr_abs" = readr::col_double(),
            "CHH_5ltr_rel" = readr::col_double(),
            "CCG_5ltr_abs" = readr::col_double(),
            "CCG_5ltr_rel" = readr::col_double(),
            "N_5ltr_abs" = readr::col_double(),
            "cn_3ltr" = readr::col_double(),
            "cn_5ltr" = readr::col_double()
        ), delim = "\t")
    } else {
        pred <- readr::read_delim(data.sheet, col_types = readr::cols(
            "ID" = readr::col_character(),
            "dfam_target_name" = readr::col_character(),
            "ltr_similarity" = readr::col_double(),
            "similarity" = readr::col_character(),
            "protein_domain" = readr::col_character(),
            "orfs" = readr::col_integer(),
            "chromosome" = readr::col_character(),
            "start" = readr::col_integer(),
            "end" = readr::col_integer(),
            "strand" = readr::col_character(),
            "width" = readr::col_integer(),
            "annotation" = readr::col_character(),
            "pred_tool" = readr::col_character(),
            "frame" = readr::col_character(),
            "score" = readr::col_character(),
            "lLTR_start" = readr::col_integer(),
            "lLTR_end" = readr::col_integer(),
            "lLTR_length" = readr::col_integer(),
            "rLTR_start" = readr::col_integer(),
            "rLTR_end" = readr::col_integer(),
            "rLTR_length" = readr::col_integer(),
            "lTSD_start" = readr::col_integer(),
            "lTSD_end" = readr::col_integer(),
            "lTSD_motif" = readr::col_character(),
            "rTSD_start" = readr::col_integer(),
            "rTSD_end" = readr::col_integer(),
            "rTSD_motif" = readr::col_character(),
            "PPT_start" = readr::col_integer(),
            "PPT_end" = readr::col_integer(),
            "PPT_motif" = readr::col_character(),
            "PPT_strand" = readr::col_character(),
            "PPT_offset" = readr::col_integer(),
            "PBS_start" = readr::col_integer(),
            "PBS_end" = readr::col_integer(),
            "PBS_strand" = readr::col_character(),
            "tRNA" = readr::col_character(),
            "tRNA_motif" = readr::col_character(),
            "PBS_offset" = readr::col_integer(),
            "tRNA_offset" = readr::col_integer(),
            "PBS/tRNA_edist" = readr::col_integer(),
            "orf.id" = readr::col_character(),
            "repeat_region_length" = readr::col_integer(),
            "PPT_length" = readr::col_integer(),
            "PBS_length" = readr::col_integer(),
            "dfam_acc" = readr::col_character(),
            "dfam_bits" = readr::col_double(),
            "dfam_e_value" = readr::col_double(),
            "dfam_bias" = readr::col_double(),
            "dfam_hmm-st" = readr::col_double(),
            "dfam_hmm-en" = readr::col_double(),
            "dfam_strand" = readr::col_character(),
            "dfam_ali-st" = readr::col_double(),
            "dfam_ali-en" = readr::col_double(),
            "dfam_env-st" = readr::col_double(),
            "dfam_env-en" = readr::col_double(),
            "dfam_modlen" = readr::col_double(),
            "dfam_target_description" = readr::col_character(),
            "Clust_Cluster" = readr::col_character(),
            "Clust_Target" = readr::col_character(),
            "Clust_Perc_Ident" = readr::col_double(),
            "Clust_cn" = readr::col_integer(),
            "TE_CG_abs" = readr::col_double(),
            "TE_CG_rel" = readr::col_double(),
            "TE_CHG_abs" = readr::col_double(),
            "TE_CHG_rel" = readr::col_double(),
            "TE_CHH_abs" = readr::col_double(),
            "TE_CHH_rel" = readr::col_double(),
            "TE_CCG_abs" = readr::col_double(),
            "TE_CCG_rel" = readr::col_double(),
            "TE_N_abs" = readr::col_double(),
            "CG_3ltr_abs" = readr::col_double(),
            "CG_3ltr_rel" = readr::col_double(),
            "CHG_3ltr_abs" = readr::col_double(),
            "CHG_3ltr_rel" = readr::col_double(),
            "CHH_3ltr_abs" = readr::col_double(),
            "CHH_3ltr_rel" = readr::col_double(),
            "CCG_3ltr_abs" = readr::col_double(),
            "CCG_3ltr_rel" = readr::col_double(),
            "N_3ltr_abs" = readr::col_double(),
            "CG_5ltr_abs" = readr::col_double(),
            "CG_5ltr_rel" = readr::col_double(),
            "CHG_5ltr_abs" = readr::col_double(),
            "CHG_5ltr_rel" = readr::col_double(),
            "CHH_5ltr_abs" = readr::col_double(),
            "CHH_5ltr_rel" = readr::col_double(),
            "CCG_5ltr_abs" = readr::col_double(),
            "CCG_5ltr_rel" = readr::col_double(),
            "N_5ltr_abs" = readr::col_double(),
            "cn_3ltr" = readr::col_double(),
            "cn_5ltr" = readr::col_double()
        ), delim = "\t")
    }
    
    
    
    return(pred)
}


