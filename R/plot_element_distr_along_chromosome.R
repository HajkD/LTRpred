#' @title Plot positions of predicted retrotransposons along chromosomes
#' @description The positionas of LTR retrotransposons predicted with \code{\link{LTRpred}}
#' will be visualized along the chromsome.
#' @param pred LTRpred.tbl generated with \code{\link{LTRpred}}.
#' @param genome.file a file path to the genome assembly file in fasta format for which chromosomes shall be visualized.
#' @param centromere_start a numeric vector storing the centromere start coordinates in the \code{genome.file}. The position in the numeric vector should correspond to the chromosome name in the \code{genome.file} fasta file. If \code{centromere_start = NULL} (default), then no centromeres will be drawn.
#'@param centromere_end a numeric vector storing the centromere end coordinates in the \code{genome.file}. The position in the numeric vector should correspond to the chromosome name in the \code{genome.file} fasta file. If \code{centromere_end = NULL} (default), then no centromeres will be drawn.
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

plot_element_distr_along_chromosome <-
  function(pred,
           genome.file,
           centromere_start = NULL,
           centromere_end = NULL,
           ...) {
    if (!file.exists(genome.file))
    stop("The genome.file '", genome.file, "' does not seem to exist. Please provide a valid path to a genome file in fasta format.", call. = FALSE)
  
  if (!is.null(centromere_start) & !is.null(centromere_end)) {
    if (!is.numeric(centromere_start))
      stop("Please provide numeric values for centromere starts.", call. = FALSE)
    if (!is.numeric(centromere_end))
      stop("Please provide numeric values for centromere end", call. = FALSE)
    
    if (length(centromere_start) != length(centromere_end))
      stop("Your centromere_start and centromere_end vectors have different lengthsm please make sure that the start and end coordinates match.", call. = FALSE)
    
    genome <- Biostrings::readDNAStringSet(genome.file)
    
    news_names <- unlist(sapply(names(genome), function(x) unlist(stringr::str_split(x," "))[1]))
    
    names(news_names) <- news_names
    
    if (length(genome) != length(centromere_start))
      stop("The number of chroosomes (= ", length(genome), ") stored in your genome file doesn't match the number of centromere start coordinates (= ",length(centromere_start), ")." , call. = FALSE)
    
  }
    
  if (is.null(centromere_start) & is.null(centromere_end)) {
      genome <- Biostrings::readDNAStringSet(genome.file)
  }
    
  if (!all(genome@ranges@NAMES == unique(pred$chromosome)))
    stop("The chromosome names from the genome.file (", paste0(genome@ranges@NAMES , collapse = ", "), ") and the chromosome names from the pred file (",paste0(unique(pred$chromosome), collapse = ", "), ") do not match. Please fix this for this function to work.", call. = FALSE)
  
  dn <- LTRpred::pred2GRanges(pred)
  GenomeInfoDb::seqlengths(dn) <- genome@ranges@width
  dn <- GenomeInfoDb::keepSeqlevels(dn, genome@ranges@NAMES)
  
  p <- ggbio::autoplot(dn, layout = "karyogram", ...)
  
  return(p)
  
  #ggbio::plotKaryogram(dn, ...) 
  
}