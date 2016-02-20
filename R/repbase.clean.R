#' @title Clean the initial Repbase database for BLAST
#' @description Clean the headers of the Repbase fasta files
#' so that headers can be used to create a blast-able database.
#' @param repbase.file fasta file storing the corresponding Repbase annotation, e.g. \code{athrep.ref}.
#' @param output.file name/path of the cleaned Repbase annotation file. 
#' @author Hajk-Georg Drost
#' @details The Repbase database can be downloaded after registration at http://www.girinst.org/repbase/.
#' The corresponding files as they are however, cannot be converted into a blast-able database.
#' Hence, a pre-filtering step is neccessary to be able to use this database with the e.g. \code{\link{repbase.query}}
#' function.  
#' @examples 
#' \dontrun{
#' # PreProcess Repbase: A thaliana
#' # and save the output into the file "Athaliana_repbase.ref"
#' repbase.clean(repbase.file = "athrep.ref",
#'               output.file  = "Athaliana_repbase.ref")
#' 
#' }
#' @references http://www.girinst.org/repbase/
#' @export
repbase.clean <- function(repbase.file, output.file){
    
    repbase <- Biostrings::readDNAStringSet(repbase.file)
    repbase@ranges@NAMES <- stringr::str_replace_all(repbase@ranges@NAMES," ","_")
    repbase@ranges@NAMES <- unlist(lapply(stringr::str_split(repbase@ranges@NAMES, "\t"), paste0, collapse = "_"))
    Biostrings::writeXStringSet(repbase,output.file)
    
}

