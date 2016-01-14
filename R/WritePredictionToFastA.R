#' @title Save the sequence of the predicted LTR Transposons in a fasta file
#' @description This function allows users to save the sequence of the predicted LTR Transposons or LTRs in a fasta file.
#' @param ltr.digest.prediction a \code{data.frame} returned by \code{\link{LTRdigest}}.
#' @param prediction.fasta.file the fasta file storing either the full LTR Transposon sequence or only the LTR sequence
#' as returned by the \code{\link{LTRdigest}} function.
#' @param output the fasta file to which the output sequences shall be stored in.
#' @author Hajk-Georg Drost
#' @details 
#' The output \code{data.frame}s returned by \code{\link{LTRdigest}} contain all information of the predicted
#' LTR retrotransposons that can be used for post-filtering steps. After these post-filtering steps
#' sequences of the remaining (filtered) candidates can be retrieved by this function.
#' @examples 
#' \dontrun{
#' # Hypothetical Example
#' # Generate LTR transposon prediction for A. thaliana
#'  LTRharvest("Genome/TAIR10_chr_all.fas")
#' 
#'  LTRdigest(input.gff3        = "TAIR10_chr_all/TAIR10_chr_all_Prediction.gff", 
#'            genome.file       = "Genome/TAIR10_chr_all.fas",
#'            trnas             = "araTha1-tRNAs.fa",
#'            hmms              = "hmm_*",
#'            cores             = 1)
#'
#' # Read the output of LTRdigest()
#' LTRdigest <- read.prediction(gff.file    = "TAIR10_chr_all_LTRdigestPrediction.gff",
#'                              tabout.file = "TAIR10_chr_all-ltrdigest_tabout.csv",
#'                              program     = "LTRdigest")
#' 
#' # Filter for LTR transposons having 100 pec. sequence similarity between their LTRs
#' FilteredLTRTransposons <- dplyr::filter(Ath.LTRdigest.prediction$ltr.retrotransposon, 
#'                                         ltr_similarity == 100)
#' 
#' # Write the sequences of these filtered LTR transposons to a fasta file
#' WritePredictionToFastA(ltr.digest.prediction = FilteredLTRTransposons, 
#'                        prediction.fasta.file = "TAIR10_chr_all-ltrdigest_complete.fas", 
#'                        output                = "AthalianaPutativeLTRTransposons.fa")
#' }
#' @seealso \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, \code{\link{read.prediction}} 
#' @export 

WritePredictionToFastA <- function(ltr.digest.prediction, prediction.fasta.file, output = "output.fa"){
    
    pred.names <- paste0(ltr.digest.prediction$sequence,"_",ltr.digest.prediction$start,"_",ltr.digest.prediction$end)
    PutativeLTRSeqs <- Biostrings::readDNAStringSet(prediction.fasta.file)
    Biostrings::writeXStringSet(PutativeLTRSeqs[match(pred.names,PutativeLTRSeqs@ranges@NAMES)],output)
    
}




