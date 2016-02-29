#' @title Perform Meta-Analyses with LTRpred
#' @description Run \code{\link{LTRpred}} on several genomes (sequentially) that are stored in a given folder.
#' @param genome.folder path to the folder storing all unmasked genomes for which \code{\link{LTRpred}}
#' based \code{de novo} LTR retrotransposon prediction shall be performed.
#' @param result.folder folder in which \code{LTRpred} results shall be stored.
#' @param ... all parameters of \code{\link{LTRpred}}.
#' @author Hajk-Georg Drost
#' @details 
#' This function provides a crawler to run \code{\link{LTRpred}} sequencially
#' on any number of genomes stored within the same folder.
#' The result will be saved on in the \code{result.folder} and can be analysed with the
#' \code{\link{meta.apply}} function.
#' @examples 
#' \dontrun{
#' # perform a meta analysis on multiple genomes
#' # stored in a genomes folder
#' LTRpred.meta(genome.folder = "Genomes/",
#'              result.folder = "LTRpredResults"
#'              trnas         = "plantRNA_Arabidopsis.fsa",
#'              hmms          = "hmm_*")
#' }
#' @export
   
LTRpred.meta <- function(genome.folder, result.folder = NULL, ...){
  
  cat("\n")
  cat("Starting LTRpred meta analysis on the following genomes: ")
  genomes <- list.files(genome.folder)
  cat("\n")
  cat("\n")
  cat(paste(genomes, sep = ", "))
  cat("\n")
  cat("\n")
  genome.names.chopped <- sapply(genomes, function(x) unlist(stringr::str_split(x, "[.]"))[1])
  
  # run meta analysis for all species sequencially
  for (i in 1:length(genomes)){
    LTRpred(...)
  }
  
  if (!is.null(result.folder)){
    # store results in result folder -> default: working directory
    file.move(from = paste0(genome.names.chopped,"_ltrpred"), 
              to   = file.path(result.folder,paste0(genome.names.chopped,"_ltrpred")))
  }
  
  cat("Finished meta analysis!")
}








