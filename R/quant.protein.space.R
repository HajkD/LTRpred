#' @title Quantify the protein space within a genome
#' @description Quantification of the protein space (= total length of
#' all annotated protein coding genes) within a given genome of interest.
#' @param file file path to a fasta file storing the protein sequences.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{quant.cds.space}}, \code{\link{quant.repeat.space}},
#' \code{\link{quant.meta.space}}
#' @export
quant.protein.space <- function(file) {
    prot_seqs <- biomartr::read_proteome(file = file)
    prot_space_in_Mbp <- sum(prot_seqs@ranges@width * 3L) / 1000000L
    nProteins <- length(prot_seqs)
    res <-
        tibble::tibble(ProtSpaceMbp = prot_space_in_Mbp, nProts = nProteins)
    return(res)
}

