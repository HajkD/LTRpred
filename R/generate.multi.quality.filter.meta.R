#' @title Run \code{quality.filter.meta} for several different ltr similarity thresholds.
#' @description A helper function to apply the \code{\link{quality.filter}} function to diverse \code{\link{LTRpred}} annotations while probing different ltr similarity thresholds.
#' @param kingdom the taxonomic kingdom of the species for which \code{\link{LTRpred}} annotations are 
#' stored in the \code{genome.folder}.
#' @param genome.folder  a file path to a folder storing the genome assembly files in fasta format that
#' were used to generate \code{\link{LTRpred}} annotations of diverse species from the same taxonomic kingdom.
#' @param ltrpred.meta.folder a file path to a folder storing \code{\link{LTRpred}} annotations of diverse species from the same taxonomic kingdom.
#' @param sim.options a numeric vector storing the ltr similarity thresholds that shall be probed.
#' @param cut.range.options a numeric vector storing the similarity cut range thresholds that shall be probed.
#' @param n.orfs minimum number of open reading frames a predicted retroelement shall possess.
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
#' @seealso \code{\link{quality.filter}}, \code{\link{quality.filter.meta}}, \code{\link{LTRpred}}, \code{\link{LTRpred.meta}}, \code{\link{read.ltrpred}}
#' @export

generate.multi.quality.filter.meta <- function(kingdom,
                                genome.folder,
                                ltrpred.meta.folder,
                                sim.options,
                                cut.range.options,
                                n.orfs = 0,
                                strategy = "default",
                                update = FALSE) {
  
  if (length(sim.options) != length(cut.range.options))
    stop("Please provide the same number of 'sim.options' and 'cut.range.options'. The length of these vectors differs ...", call. = FALSE)
  
  ltrpred_meta_result_i <- vector("list", length(sim.options))
    
  for (i in seq_len(length(sim.options)))  {
    
    ltrpred_meta_result_i[i] <- list(quality.filter.meta(
      kingdom = paste0(kingdom,"_", sim.options[i]),
      genome.folder = genome.folder,
      ltrpred.meta.folder = ltrpred.meta.folder,
      sim = sim.options[i],
      cut.range = cut.range.options[i],
      n.orfs = n.orfs,
      strategy = strategy,
      update = update
    ))
  }
  
  names(ltrpred_meta_result_i) <- paste0(kingdom,"_", sim.options)
  return(ltrpred_meta_result_i)
}