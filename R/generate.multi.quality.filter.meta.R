#' @title Run \code{quality.filter.meta} for several different ltr similarity thresholds.
#' @description 
#' @param kingdom
#' @param genome.folder
#' @param ltrpred.meta.folder
#' @param sim.options
#' @param cut.range.options
#' @param n.orfs
#' @param strategy
#' @param update
#' @author Hajk-Georg Drost
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