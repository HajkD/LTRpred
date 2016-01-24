#' @title Read output of USEARCH's findorfs output
#' @description This function reads the output of the \code{USEARCH}
#' \code{fastx_findorfs} command line tool and stores the 
#' sequence id and number of predicted ORFs in a \code{\link{data.frame}} object.
#' @param findorfs.file
#' @author Hajk-Georg Drost
#' @export
read.findorfs.output <- function(findorfs.file){
  
  ReadSeqFile <- Biostrings::readDNAStringSet(findorfs.file)
  SeqFile.table <- table(sapply(ReadSeqFile@ranges@NAMES, 
                                function(x) unlist(stringr::str_split(x, "[|]"))[1]))
  ORFCount.df <- dplyr::data_frame(seq.id = names(SeqFile.table), 
                                   orfs = as.numeric(SeqFile.table))
  GenomicLocus <- as.data.frame(do.call(rbind, sapply(ORFCount.df$seq.id, function(x){
    as.numeric(unlist(stringr::str_split(unlist(stringr::str_split(x,"__"))[2],"_")))
  })), row.names = FALSE)
  names(GenomicLocus) <- c("start","end")
  remove.NA <- which(is.na(GenomicLocus$start) | is.na(GenomicLocus$end))
  GenomicLocus <- GenomicLocus[-remove.NA, ]
  ORFCount.df <- ORFCount.df[-remove.NA, ]
  ORFCount.df <- dplyr::mutate(ORFCount.df, start = unlist(GenomicLocus$start), end = unlist(GenomicLocus$end))
  
  return (ORFCount.df) 
}






