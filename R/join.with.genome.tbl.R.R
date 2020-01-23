#' @title Join \code{gm_files} returned by \code{generate.multi.quality.filter.meta} with a genome information table
#' @description Join \code{gm_files} returned by \code{generate.multi.quality.filter.meta} with a genome information table. 
#' @param meta.summary.file a meta.summary.file.
#' @param genome.tbl a genome.tbl.
#' @author Hajk-Georg Drost
#' @export

join.with.genome.tbl <- function(meta.summary.file, genome.tbl) {
  
  for (i in seq_len(length(meta.summary.file))) {
    meta.summary.file[[i]]$gm_file <-
      dplyr::left_join(meta.summary.file[[i]]$gm_file, genome.tbl, by = "organism")
    meta.summary.file[[i]]$gm_file <- dplyr::filter(meta.summary.file[[i]]$gm_file, !is.na(order))
  }
  
  return(meta.summary.file)
}