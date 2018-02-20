#' @title Quantify the total repeat space from Repeat Masker output in Mbp
#' @description 
#' @param file path to Repeat Masker output file.
#' @param repeat.type type of element for which total repeat space shall be quantified.
#' Options are:
#' 
#' \itemize{
#'  \item \code{repeat.type = "all"} include all types of repeats.
#'  \item \code{repeat.type = "LTR"} include only LTR retrotransposons.
#' }
#' @seealso \code{\link{meta.seq.space}}
#' @author Hajk-Georg Drost
#' @export

totalRepeatSpace <- function(file, repeat.type = "LTR") {
  
  if (!is.element(repeat.type, c("LTR", "all")))
    stop("Please specify a type argument that is supported by this function.", call. = FALSE)
  
  if (!file.exists(file))
    stop("The file '",file,"' could not be found. Please provide a valid path to a Repeat Masker output file.", call. = FALSE)
  
  message("\n")
  message("Importing Repeat Masker Annotation from '",file,"' ...")
  repeat_file <- biomartr::read_rm(file = file)
  
  if (repeat.type == "LTR") {
    message("Selecting only LTR retrotransposon associated repeats ...")
    matching_class <- NULL
    repeat_file <- dplyr::filter(repeat_file, stringr::str_detect(matching_class, "LTR"))
    
    message("Computing the total LTR retrotransposon space in Mbp ...")
    repeat_space_in_Mbp <- sum(repeat_file$qry_width, na.rm = TRUE) / 1000000L
    names(repeat_space_in_Mbp) <- "ltr_retro_space_mbp"
  }
  
  if (repeat.type == "all") {
    message("Computing the total repeat space in Mbp ...")
    repeat_space_in_Mbp <- sum(repeat_file$qry_width, na.rm = TRUE) / 1000000L
    names(repeat_space_in_Mbp) <- "total_repeat_space_mbp"
  }
  
  message("Finished!")
  return(repeat_space_in_Mbp)
}
