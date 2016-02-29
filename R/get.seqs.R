#' @title Quickly retrieve the sequences of a \code{Biostrings} object
#' @description Helper function to retrieve the sequences of a \code{Biostrings} object.
#' @param object a \code{Biostrings} object.
#' @param type either \code{type = "DNA"} or \code{type = "AA"}.
#' @author Hajk-Georg Drost
#' @examples 
#' # read example sequences
#' seqs <- system.file("nt.fa",package = "LTRpred")
#' input.seqs <- read.seqs(seqs)
#' 
#' retrieve sequences of Biostrings object
#' head(get.seqs(input.seqs))
#' @export

get.seqs <- function(object, type = "DNA"){
  
  if (!is.element(type,c("AA","DNA")))
    stop ("Please select a valid type: either DNA or AA.")
  
  if (type == "DNA")
    return (lapply(object, function(x) Biostrings::DNAString(x)))
  if (type == "AA")
    return (lapply(object, function(x) Biostrings::AAString(x)))
}


