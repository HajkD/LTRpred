#' @title Low level function to detect motifs in strings
#' @description Find a specific motif or a sequence of motifs within
#' genomic sequences.
#' @param seq.file path to the genomic sequecne file of interest (e.g. LTR TE seqs predicted by
#' \code{\link{LTRpred}}).
#' @param motif a character string or vector of strings which shall be counted within each sequence.
#' @param as.ratio shall count values be returned as asbolute frequency (count value) or as relative frequency (percentage).
#' @author Hajk-Georg Drost
#' @examples 
#' # find number of "CG" motifs in predicted LTR transposons
#' motif.count(seq.file = system.file("LTRseqs.fas",package = "LTRpred"), 
#'             motif    = "CG")
#'             
#' # find number of "CG" motifs in predicted LTR transposons: rel. frequency
#' motif.count(seq.file = system.file("LTRseqs.fas",package = "LTRpred"), 
#'             motif    = "CG",
#'             as.ratio = TRUE)             
#' @export
       
motif.count <- function(seq.file, motif, as.ratio = FALSE){
   
  # read sequence
  seqs <- Biostrings::readDNAStringSet(seq.file)
  
  if (as.ratio){
    res <- Biostrings::vcountPattern(motif,seqs) / Biostrings::nchar(seqs)
    names(res) <- seqs@ranges@NAMES
  }
   
  if (!as.ratio){
    res <- Biostrings::vcountPattern(motif,seqs)
    names(res) <- seqs@ranges@NAMES
  }
   
  return (res)
}






