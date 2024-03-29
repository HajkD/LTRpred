#' @title Format LTR prediction data to \code{GRages} object
#' @description This function formats the LTR prediction \code{\link{data.frame}}
#' generated by \code{\link{LTRpred}} to a \code{GRages} object.
#' @param LTR.data the LTR prediction \code{\link{data.frame}} generated by \code{\link{LTRpred}}.
#' @param similarity.threshold the LTR similarity threshold that shall be used to define young and old retrotransposons.
#' @author Hajk-Georg Drost
#' @details 
#' The \code{GRages} object is defined by the following columns:
#' \itemize{
#' \item \code{chromosome}
#' \item \code{start}
#' \item \code{end}
#' \item \code{strand}
#' \item \code{element_name}
#' \item \code{ltr_similarity}
#' \item \code{orfs}
#' \item \code{age}
#' }
#' @export
pred2GRanges <- function(LTR.data, similarity.threshold = 95){
  
  res <- dplyr::data_frame(
                           chromosome = LTR.data$chromosome, 
                           start      = LTR.data$start, 
                           end        = LTR.data$end,
                           strand     = LTR.data$strand,
                           element_name = LTR.data$ID,
                           ltr_similarity = LTR.data$similarity,
                           orfs           = LTR.data$orfs,
                           age            = ifelse(LTR.data$ltr_similarity  < similarity.threshold, "old", "young"))
  
  res_gr <- GenomicRanges::makeGRangesFromDataFrame(res, keep.extra.columns = TRUE)
  
  return(res_gr)
}


