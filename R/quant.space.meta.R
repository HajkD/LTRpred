#' @title Quantify the genomic loci space for multiple genomes, proteomes or Repeat Masker annotation files
#' @description Quantification of the genomic loci space (= total length of
#' all annotated genomic features) within a given genome of interest.
#' @param folder path to the folder storing the genomes, proteomes, or 
#' Repeat Masker files of interest.
#' @param type The following genomic features can be quantified:
#' \itemize{
#' \item \code{type = 'cds'}: coding sequence files in fasta format (see \code{\link[biomartr]{read_genome}}).
#' \item \code{type = 'proteome'}: proteome files in fasta format (see \code{\link[biomartr]{read_proteome}}).
#' \item \code{type = 'rm'}: Repeat Masker annotation files in fasta format (see \code{\link[biomartr]{read_rm}}).
#' }
#' @author Hajk-Georg Drost
#' @export
quant.space.meta <- function(folder, type = "proteome") {
    if (!is.element(type, c("proteome", "cds", "rm")))
        stop(
            "Please select a sequence type for which a sequence space can",
            " be computed, e.g. 'cds', 'proteome' or 'rm'."
            
        )
    
    if (type == "proteome") {
        # retrieve sequence file names and remove doc_* files
        get_files <- list.files(folder)
        get_files <-
            get_files[!unlist(stringr::str_detect(get_files, "doc_"))]
        get_files <- file.path(folder, get_files)
        
        prot_space_all <- dplyr::bind_rows(lapply(get_files, function(x) {
            message("Processing file '",x,"' ...")
            prot_space <- quant.protein.space(file = x)
            return(prot_space)
        }))
        
        res <-
            tibble::tibble(organism = basename(get_files))
        res <- dplyr::bind_cols(res, prot_space_all)
        
        return(res)
    }
    
    
    if (type == "rm") {
        # retrieve sequence file names and remove doc_* files
        get_files <- list.files(folder)
        get_files <-
            get_files[!unlist(stringr::str_detect(get_files, "doc_"))]
        get_files <- file.path(folder, get_files)
        
        repeat_space_all <- dplyr::bind_rows(lapply(get_files, function(x) {
            message("Processing file '",x,"' ...")
            repeat_space <- quant.repeat.space(file = x)
            return(repeat_space)
        }))
        
        res <-
            tibble::tibble(organism = basename(get_files))
        res <- dplyr::bind_cols(res, repeat_space_all)
        
        return(res)
    }
    
    if (type == "cds") {
        # retrieve sequence file names and remove doc_* files
        get_files <- list.files(folder)
        get_files <-
            get_files[!unlist(stringr::str_detect(get_files, "doc_"))]
        get_files <- file.path(folder, get_files)
        
        cds_space_all <- dplyr::bind_rows(lapply(get_files, function(x) {
            message("Processing file '",x,"' ...")
            cds_space <- quant.cds.space(file = x)
            return(cds_space)
        }))
        
        res <-
            tibble::tibble(organism = basename(get_files))
        res <- dplyr::bind_cols(res, cds_space_all)
        
        return(res)
    }
}
