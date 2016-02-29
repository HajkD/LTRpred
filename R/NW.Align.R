# Needleman-Wunsch Global Alignment function
NW.Align <- function(seq1,seq2, ...){
  return (Biostrings::pairwiseAlignment(pattern = seq1, subject = seq2,type = "global", ...))
}