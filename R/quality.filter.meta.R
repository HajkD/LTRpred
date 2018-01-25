#' @title Pipeline to eliminate false positive predictions of retrotransposons at a metagenomic scale
#' @description This function takes the file paths to the genomes folder and \code{\link{LTRpred.meta}} output folder as input and eliminates false positive retrotransposon predictions on a metagenomic scale.
#' @param kingdom a character string specifying the kingdom of life to which genomes annotated with \code{\link{LTRpred.meta}} belong to. E.g. \code{kingdom = "Plants"}. If annotates of a variety of kingdoms have been done users can for example specify \code{kingdom = "Various"}.
#' @param genome.folder path to folder storing the genome assembly files that were used for \code{\link{LTRpred.meta}} predictions.
#' @param ltrpred.meta.folder path to folder storing the \code{\link{LTRpred.meta}} output files.
#' @param sim LTR similarity threshold. Only putative LTR transposons that fulfill this 
#' LTR similarity threshold will be retained.
#' @param cut.range a numeric number indicating the interval size for binning LTR similarities.
#' @param n.orfs minimum number of ORFs detected in the putative LTR transposon.
#' @param strategy quality filter strategy. Options are
#' \itemize{
#' \item \code{strategy = "default"} : see section \code{Quality Control} 
#' \item \code{strategy = "stringent"} : in addition to filter criteria specified in section \code{Quality Control},
#' the filter criteria \code{!is.na(protein_domain)) | (dfam_target_name != "unknown")} is applied
#' }
#' @param update shall already existing \code{_SimilarityMatrix.csv} and \code{_GenomeInfo.csv} files be updated (\code{update = TRUE}) or can the already existing files be used (\code{update = FALSE})?
#' @author Hajk-Georg Drost
#' @details 
#' \strong{Quality Control}
#' 
#' \itemize{
#' \item \code{ltr.similarity}: Minimum similarity between LTRs. All TEs not matching this
#'  criteria are discarded.
#'  \item \code{n.orfs}: minimum number of Open Reading Frames that must be found between the
#'   LTRs. All TEs not matching this criteria are discarded.
#'  \item \code{PBS or Protein Match}: elements must either have a predicted Primer Binding
#'  Site or a protein match of at least one protein (Gag, Pol, Rve, ...) between their LTRs. All TEs not matching this criteria are discarded.
#'  \item The relative number of N's (= nucleotide not known) in TE <= 0.1. The relative number of N's is computed as follows: absolute number of N's in TE / width of TE.
#' }
#' @return A list with to list elements \code{sim_file} and \code{gm_file}. Each list element stores a \code{data.frame}:
#'   \itemize{
#'   \item \code{sim_file} (similarity file)
#'          \itemize{
#'          This \code{data.frame} stores the information 
#'                   \item 
#'                   }
#'    \item \code{gm_file} (genome metrics file)
#'          \itemize{
#'          This \code{data.frame} stores the information
#'                   \item
#'                   }
#'    }               
#' @seealso \code{\link{LTRpred}}, \code{\link{LTRpred.meta}}, \code{\link{read.ltrpred}} 
#' @export 

quality.filter.meta <-
  function(kingdom,
           genome.folder,
           ltrpred.meta.folder,
           sim,
           cut.range = 2,
           n.orfs,
           strategy,
           update = FALSE) {
    
    
    sim_file <- paste0(kingdom, "_SimilarityMatrix.csv")
    gm_file <- paste0(kingdom, "_GenomeInfo.csv")
    
    if (!(file.exists(sim_file) & file.exists(gm_file) & (!update))) {
      genome.summary(
        genome.folder       = genome.folder, # path to folder storing genome assemblies of species
        ltrpred.meta.folder = ltrpred.meta.folder, # path to result folder 
        quality.filter      = TRUE, # apply filter to reduce false positives ->
        sim                 = sim, # LTR retrotransposons should have >= 70% seq. similarity between their LTRs
        cut.range           = cut.range,
        n.orfs              = n.orfs, # LTR retrotransposons should have >= 0 Open Reading Frame
        strategy            = strategy,
        file.name           = kingdom
      )
    } else {
      message("\n")
      message("The files '",sim_file,"' and '",gm_file,"' were already found in the current working directory and will be used for further processing. If you wish to re-run the quality filtering and generate new '",sim_file,"' and '",gm_file,"' files please specify the argument 'update = TRUE'.")
    }
  
    suppressMessages(sim_file_import <-
      readr::read_delim(sim_file, delim = ";"))
    
    gm_file_import <-
      readr::read_delim(gm_file, delim = ";", col_types = readr::cols(
          organism = readr::col_character(),
          n_ltrs = readr::col_integer(),
          total_ltrs_nucl_mbp = readr::col_double(),
          total_ltrs_nucl_freq = readr::col_double(),
          n_ltrs_freq = readr::col_double(),
          genome_size_nucl_mbp = readr::col_double(),
          NNN_freq = readr::col_double()
      ))
    
    sim_file_import <- dplyr::mutate(sim_file_import, kingdom = rep(kingdom, nrow(sim_file_import)))
    sim_file_import <- dplyr::select(sim_file_import, kingdom, 1:(ncol(sim_file_import) - 1))
    gm_file_import <- dplyr::mutate(gm_file_import, kingdom = rep(kingdom, nrow(gm_file_import)))
    gm_file_import <- dplyr::select(gm_file_import, kingdom, 1:(ncol(gm_file_import) - 1))
    
    
    res <- list(sim_file = sim_file_import, gm_file = gm_file_import) 
    
    return(res)
   }




