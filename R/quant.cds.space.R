#' @title Quantify the coding sequence space within a genome
#' @description Quantification of the cds space (= total length of
#' all annotated coding sequences) within a given genome of interest.
#' @param file file path to a fasta file storing the cds sequences.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{quant.protein.space}}, \code{\link{quant.repeat.space}}
#' @export
quant.cds.space <- function(file) {
    cds_seqs <- biomartr::read_cds(file = file)
    cds_space_in_Mbp <- sum(cds_seqs@ranges@width) / 1000000L
    nCDSs <- length(cds_seqs)
    res <-
        tibble::tibble(CDSSpaceMbp = cds_space_in_Mbp, nCDSs = nCDSs)
    return(res)
}

