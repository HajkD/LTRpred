#' @title Helper function to add species names to headers within fasta files
#' @description Reads a fasta file and adds the defined \code{species} name to each
#' header entry of the input fasta file.
#' @param file input fasta file.
#' @param species a character string specifying the species name to be added to the headers. 
#' Format will be: \code{species_}*, where * stands for the original header.
#' @param output a character string denoting the name of the renamed output fasta file.
#' @param append logical value. If \code{TRUE} output will be appended to file; otherwise, it will overwrite the contents of file.
#' @author Hajk-Georg Drost
#' @return Writes a new fasta file with renamed headers.
#' @export

rename.fasta <- function(file, species, output = "renamed_fasta.fa", append = FALSE){
  
  fa.file <- Biostrings::readDNAStringSet(file)
  fa.file@ranges@NAMES <-
    unlist(sapply(fa.file@ranges@NAMES, function(x) {
      if (stringr::str_detect(x, "\\["))
        x <- stringr::str_replace(x, "\\[", "")
      
      stringr::str_replace(x, x, paste0(species, "_", x))
    }))
  Biostrings::writeXStringSet(fa.file, output, append = append)
}
