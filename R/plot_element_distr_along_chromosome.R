#' @title Plot positions of predicted retrotransposons along chromosomes
#' @description The positionas of LTR retrotransposons predicted with \code{\link{LTRpred}}
#' will be visualized along the chromsome.
#' @param pred LTRpred.tbl generated with \code{\link{LTRpred}}.
#' @param genome.file a file path to the genome assembly file in fasta format for which chromosomes shall be visualized.
#' @param centromere_starts a numeric vector containing the start coordinates of the centromers (for all chromosomes in \code{genome.file}).
#' @param ... additional arguments that shall be passed to the visualization function \code{\link[ggbio]{autoplot}}.
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' test_genome <- system.file("Hsapiens_ChrY.fa", package = "LTRpred")
#' test_pred <- LTRpred::read.ltrpred(
#' system.file("Hsapiens_ChrY_LTRpred_DataSheet.tsv", 
#'              package = "LTRpred"))
#' test_centromere_starts <- 55000
#' # generate visualization
#' LTRpred::plot_element_distr_along_chromosome(test_pred, test_genome, test_centromere_starts)
#' }
#' @export

plot_element_distr_along_chromosome <- function(pred, genome.file, centromere_starts, ...){
  
  if (!file.exists(genome.file))
    stop("The genome.file '", genome.file, "' does not seem to exist. Please provide a valid path to a genome file in fasta format.", call. = FALSE)
  
  if (!is.numeric(centromere_starts))
    stop("Please provide numeric values for centromere starts.", call. = FALSE)
  
  genome <- Biostrings::readDNAStringSet(genome.file)
  
  if (length(genome) != length(centromere_starts))
    stop("The number of chroosomes (= ", length(genome), ") stored in your genome file doesn't match the number of centromere start coordinates (= ",length(centromere_starts), ")." , call. = FALSE)
  

  if (!all(genome@ranges@NAMES == unique(pred$chromosome)))
    stop("The chromosome names from the genome.file (", paste0(genome@ranges@NAMES , collapse = ", "), ") and the chromosome names from the pred file (",paste0(unique(pred$chromosome), collapse = ", "), ") do not match. Please fix this for this function to work.", call. = FALSE)
  
  dn <- LTRpred::pred2GRanges(pred)
  GenomeInfoDb::seqlengths(dn) <- genome@ranges@width
  dn <- GenomeInfoDb::keepSeqlevels(dn, genome@ranges@NAMES)
  
  ggbio::autoplot(dn, layout = "karyogram", ...)
  
  #ggbio::plotKaryogram(dn, ...) 
  
}