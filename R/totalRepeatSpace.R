#' @title Quantify the total repeat space from Repeat Masker output in Mbp
#' @description 
#' @param file path to Repeat Masker output file.
#' @param type type of element for which total repeat space shall be quantified.
#' Options are:
#' 
#' \itemize{
#'  \item \code{type = "all"} include all types of repeats.
#'  \item \code{type = "LTR"} include only LTR retrotransposons.
#' }
#' @author Hajk-Georg Drost
#' @export

totalRepeatSpace <- function(file, type = "LTR") {
  
  if (!is.element(type, c("LTR", "all")))
    stop("Please specify a type argument that is supported by this function.", call. = FALSE)
  
  if (!file.exists(file))
    stop("The file '",file,"' could not be found. Please provide a valid path to a Repeat Masker output file.", call. = FALSE)
  
  message("Importing Repeat Masker Annotation from '",file,"' ...")
  repeat_file <- biomartr::read_rm(file = file)
  
  if (type == "LTR") {
    message("Selecting only LTR retrotransposon associated repeats ...")
    matching_class <- NULL
    repeat_file <- dplyr::filter(repeat_file, stringr::str_detect(matching_class, "LTR"))
  }
  
  message("Computing the total repeat space in Mbp ...")
  repeat_space_in_Mbp <- sum(repeat_file$qry_width, na.rm = TRUE) / 1000000L
  names(repeat_space_in_Mbp) <- "repeat_space_mbp"
  message("Finished!")
  return(repeat_space_in_Mbp)
}
