#' @title Open Reading Frame Prediction
#' @description This function provides a wrapper to the \code{USEARCH}
#' \code{fastx_findorfs} command line tool to predict ORFs in
#' a given input fasta file and read the output as \code{\link{data.frame}} object.
#' @param seq.file a fasta file storing the sequences for which open reading frames shall be predicted.
#' @param orf.style type of predicting open reading frames (see documentation of USEARCH).
#' @param min.codons minimum number of codons in the predicted open reading frame.
#' @param trans.seqs logical value indicating wheter or not predicted open reading frames
#' shall be translated and the corresponding protein sequences stored in the output folder. 
#' @param output path to the folder in which predicted open reading frame sequences shall be stored.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, \code{\link{LTRpred}}
#' @references 
#' Robert Edgar. Search and clustering orders of magnitude faster than BLAST. Bioinformatics (2010) 26 (19): 2460-2461.
#' @export

ORFpred <- function(seq.file, 
                    orf.style  = 7, 
                    min.codons = 200, 
                    trans.seqs = FALSE,
                    output     = NULL){
   
    test_installation_usearch()
    
    if (!trans.seqs) {
        if (is.null(output)) {
            system(
                paste0(
                    "usearch -fastx_findorfs ",
                    ws.wrap.path(seq.file),
                    " -ntout ",
                    paste0(basename(seq.file), "_ORF_prediction_nt.fsa"),
                    " -orfstyle ",
                    orf.style,
                    " -mincodons ",
                    min.codons
                )
            )
        } else {
            system(
                paste0(
                    "usearch -fastx_findorfs ",
                    ws.wrap.path(seq.file),
                    " -ntout ",
                    file.path(output, paste0(basename(seq.file), "_ORF_prediction_nt.fsa")),
                    " -orfstyle ",
                    orf.style,
                    " -mincodons ",
                    min.codons
                )
            )
        }
        
    } else {
        if (is.null(output)) {
            system(
                paste0(
                    "usearch -fastx_findorfs ",
                    ws.wrap.path(seq.file),
                    " -ntout ",
                    paste0(basename(seq.file), "_ORF_prediction_nt.fsa"),
                    " -aaout ",
                    paste0(basename(seq.file), "_ORF_prediction_aa.fsa"),
                    " -orfstyle ",
                    orf.style,
                    " -mincodons ",
                    min.codons
                )
            )
        } else {
            system(
                paste0(
                    "usearch -fastx_findorfs ",
                    ws.wrap.path(seq.file),
                    " -ntout ",
                    file.path(output, paste0(
                        basename(seq.file), "_ORF_prediction_nt.fsa"
                    )),
                    " -aaout ",
                    file.path(output, paste0(
                        basename(seq.file), "_ORF_prediction_aa.fsa"
                    )),
                    " -orfstyle ",
                    orf.style,
                    " -mincodons ",
                    min.codons
                )
            )
        }
    }
    
    if (!is.null(output)) {
        orf.file <- file.path(output, paste0(basename(seq.file), "_ORF_prediction_nt.fsa"))
        if (file.info(orf.file)$size == 0 ||
            is.na(file.info(orf.file)$size == 0)) {
            ReadSeqFile <- Biostrings::readDNAStringSet(seq.file)
            ORFCount.df <-
                dplyr::data_frame(seq.id = ReadSeqFile@ranges@NAMES,
                                  orfs = rep(0, length(ReadSeqFile@ranges@NAMES)))
        } else {
            ORFCount.df <- read.orfs(orf.file)
        }
    } else {
        orf.file <- paste0(basename(seq.file), "_ORF_prediction_nt.fsa")
        if (file.info(orf.file)$size == 0 ||
            is.na(file.info(orf.file)$size == 0)) {
            ReadSeqFile <- Biostrings::readDNAStringSet(seq.file)
            ORFCount.df <-
                dplyr::data_frame(seq.id = ReadSeqFile@ranges@NAMES,
                                  orfs = rep(0, length(ReadSeqFile@ranges@NAMES)))
        } else {
            ORFCount.df <- read.orfs(orf.file)
        }
    }
    
    return(ORFCount.df)
}

