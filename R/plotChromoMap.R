#' @title Plot a predicted retrotransposons along the chromosomes
#' @description LTR retrotransposons predicted with \code{\link{LTRpred}}
#' will be visualized along the chromsome.
#' @param pred LTRpred.tbl generated with \code{\link{LTRpred}}.
#' @param genome.file a file path to the genome assembly file in fasta format for which chromosomes shall be visualized.
#' @param centromere_starts a numeric vector containing the start coordinates of the centromers (for all chromosomes in \code{genome.file}).
#' @param \ellipsis
#' @author Hajk-Georg Drost
#' @examples 
#' test_genome <- system.file("Hsapiens_ChrY.fa", package = "LTRpred")
#' test_pred <- LTRpred::read.ltrpred(system.file("Hsapiens_ChrY_LTRpred_DataSheet.tsv", package = "LTRpred"))
#' test_centromere_starts <- 55000
#' # generate visualization
#' LTRpred::plotChromoMap(test_pred, test_genome, test_centromere_starts)
#' @references https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html
#' @export

plotChromoMap <- function(pred, genome.file, centromere_starts, ...){
  
  if (!file.exists(genome.file))
    stop("The genome.file '", genome.file, "' does not seem to exist. Please provide a valid path to a genome file in fasta format.", call. = FALSE)
  
  if (!is.numeric(centromere_starts))
    stop("Please provide numeric values for centromere starts.", call. = FALSE)
  
  genome <- Biostrings::readDNAStringSet(genome.file)
  
  if (length(genome) != length(centromere_starts))
    stop("The number of chroosomes (= ", length(genome), ") stored in your genome file doesn't match the number of centromere start coordinates (= ",length(centromere_starts), ")." , call. = FALSE)
  
  chromoMap_chromosome_data <- dplyr::data_frame( chromosome_name = genome@ranges@NAMES,
                                                  chromosome_start = genome@ranges@start,
                                                  chromosome_end = genome@ranges@width,
                                                  centromere_start = centromere_starts)
  
  if (!all(genome@ranges@NAMES == unique(pred$chromosome)))
    stop("The chromosome names from the genome.file (", paste0(genome@ranges@NAMES , collapse = ", "), ") and the chromosome names from the pred file (",paste0(unique(pred$chromosome), collapse = ", "), ") do not match. Please fix this for this function to work.", call. = FALSE)
    
  chromoMap_chromosome_file <- file.path(tempdir(), paste0("chromoMap_chromosome_file_", unlist(stringr::str_split(basename(genome.file), "[.]"))[1], ".txt"))
  chromoMap_annotation_file <- file.path(tempdir(), paste0("chromoMap_annotation_file_", unique(pred$species), ".txt"))
  
  chromoMapChromosomeFile(data = chromoMap_chromosome_data, output = chromoMap_chromosome_file)
  pred2chromoMap(LTR.data = pred, output = chromoMap_annotation_file)
  
  chromoMap::chromoMap(chromoMap_chromosome_file, chromoMap_annotation_file, ...)
}