#' @title Import sequences of LTRharvest predicted LTR transposons
#' @description This function 
#' @param ltr.fasta.file
#' @author Hajk-Georg Drost
#' @export

ReadLTRharvestPredictionSeqs <- function(ltr.fasta.file){
    
    
    PredictedLTRSeqs <- Biostrings::readDNAStringSet(input,"fasta")
    HeaderInformation <- PredictedLTRSeqs@ranges@NAMES
    SeqInformation <- do.call(rbind,sapply(HeaderInformation, function(x) noquote(stringr::str_split(stringr::str_replace(stringr::str_replace(stringr::str_extract(x,"[?<=\\[].*?[?=\\]]"),"\\[",""),"\\]",""),","))))
    colnames(SeqInformation) <- c("start","end")
    ChrID <- sapply(rownames(SeqInformation), function(y) stringr::str_split(y, " \\(")[[1]][1])
    
    SeqInformation.df <- data.frame(chromosome = ChrID, start = as.numeric(SeqInformation[ , "start"]), end = as.numeric(SeqInformation[ , "end"]))
    SeqInformation.df <- dplyr::mutate(SeqInformation.df, width = (end - start) + 1)

    return(SeqInformation.df)
}
