#' @title Combine LTRpred prediction files
#' @description Given a path to a folder that stores several LTRpred prediction files,
#' this function imports the individual prediction files and combines them to one
#' large LTRpred prediction file.
#' @param folder path to folder in which prediction files are stored.
#' @param quality.filter shall false positives be filtered out as much as possible or not.
#' @param sim If \code{quality.filter = TRUE}: LTR similarity threshold. Only putative LTR transposons that fulfill this 
#' LTR similarity threshold will be retained.
#' @param n.orfs If \code{quality.filter = TRUE}: minimum number of ORFs detected in the putative LTR transposon.
#' @author Hajk-Georg Drost
#' @export
combinePreds <- function(folder, quality.filter = FALSE, sim = 70, n.orfs = 0) {
  input.files <- list.files(folder)
  q.tibble.list <- vector("list", length(input.files))
  
  for (i in seq_along(input.files)) {
    if (quality.filter) {
      q.tibble.list[i] <-
        list(LTRpred::quality.filter(
          LTRpred::read.ltrpred(file.path(folder, input.files[i])),
          sim = sim,
          n.orfs = n.orfs
        ))
    } else {
      q.tibble.list[i] <-
        list(LTRpred::read.ltrpred(file.path(folder, input.files[i])))
    }
  }
  return(dplyr::bind_rows(q.tibble.list))
}