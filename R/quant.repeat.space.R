#' @title Quantify the repeat space within a genome
#' @description Quantification of the repeat space (= total length of
#' all Repeat Masker annotated repeats) within a given genome of interest.
#' @param file file path to a fasta file storing the Repeat Masker annotation file.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{quant.protein.space}}, \code{\link{quant.cds.space}},
#' \code{\link{quant.meta.space}}
#' @export
quant.repeat.space <- function(file) {
    repeat_file <- biomartr::read_rm(file = file)
    repeat_space_in_Mbp <-
        sum(repeat_file$qry_width, na.rm = TRUE) / 1000000L
    nRepeats <- length(repeat_file)
    res <-
        tibble::tibble(RepeatSpaceMbp = repeat_space_in_Mbp, nRepeats = nRepeats)
    return(res)
}
