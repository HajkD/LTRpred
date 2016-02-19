#' @title Import sequences of predicted LTR transposons
#' @description This function 
#' @param seq.file a sequence file in fasta format. 
#' @author Hajk-Georg Drost
#' @export

read.seqs <- function(seq.file, program = "LTRharvest"){
    
    if (!is.element(program, c("LTRharvest")))
      stop ("Please select a prediction program that is supported by this function.")
  
    PredictedLTRSeqs <- Biostrings::readDNAStringSet(seq.file,"fasta")
    HeaderInformation <- PredictedLTRSeqs@ranges@NAMES
    SeqInformation <- do.call(rbind,sapply(HeaderInformation, function(x) noquote(stringr::str_split(stringr::str_replace(stringr::str_replace(stringr::str_extract(x,"[?<=\\[].*?[?=\\]]"),"\\[",""),"\\]",""),","))))
    colnames(SeqInformation) <- c("start","end")
    ChrID <- sapply(rownames(SeqInformation), function(y) stringr::str_split(y, " \\(")[[1]][1])
    
    SeqInformation.df <- data.frame(chromosome = ChrID, start = as.numeric(SeqInformation[ , "start"]), end = as.numeric(SeqInformation[ , "end"]))
    SeqInformation.df <- dplyr::mutate(SeqInformation.df, width = (end - start) + 1)

    return(SeqInformation.df)
}
