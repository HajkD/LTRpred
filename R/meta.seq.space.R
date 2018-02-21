#' @title Quantification of the repeat space for multiple \code{Repeat Masker} files
#' @description The repeat space for multiple \code{Repeat Masker} files is quantified and
#' stored in one \code{tibble}.
#' @param rm.folder file path to the \code{Repeat Masker} files.
#' @param type type of input data. Option are:
#' \itemize{
#'  \item \code{type = "rm"} - \code{Repeat Masker} files.
#' }
#' @param repeat.type type of element for which total repeat space shall be quantified.
#' Options are:
#' 
#' \itemize{
#'  \item \code{repeat.type = "all"} include all types of repeats.
#'  \item \code{repeat.type = "LTR"} include only LTR retrotransposons.
#' }
#' @param rename logical value indicating whether or not long file names shall be transformed to short scientific names.
#' @param save.file logical value indicating whether or not the result \code{tibble} shall be stored locally.
#' @author Hajk-Georg Drost
#' @export

meta.seq.space <-
  function(rm.folder,
           type = "rm",
           repeat.type = "all",
           rename = FALSE,
           save.file = TRUE
           ) {
    
  if (!is.element(type, c("rm")))
    stop(
      "Please select a sequence type for which a sequence space can",
      " be computed, e.g. 'genome', 'proteome' or 'rm'.", call. = FALSE)
  
  if (!file.exists(rm.folder))
    stop("The folder path '",rm.folder,"' does not exists. Please provide a valid path to a folder storing Repeat Masker files.", call. = FALSE)
  
  if (type == "rm") {
    # retrieve sequence file names and remove doc_* files
    get_files <- list.files(rm.folder)
    get_files <-
      get_files[!unlist(stringr::str_detect(get_files, "doc_"))]
    get_files <- file.path(rm.folder, get_files)
    
    repeat_space_all <- unlist(lapply(get_files, function(x) {
      message("Processing file '",x,"' ...")
      repeat_space <- totalRepeatSpace(file = x, repeat.type = repeat.type)
      return(repeat_space)
    }))
    
    if (repeat.type == "all") {
      res <-
        tibble::tibble(organism = basename(get_files),
                       total_repeat_space_mbp = repeat_space_all)
    }
    
    if (repeat.type == "ltr") {
      res <-
        tibble::tibble(organism = basename(get_files),
                       ltr_retro_space_mbp = repeat_space_all)
    }
    
    if (rename)
      res <- assign.short.sci.name(res)
    
    if (save.file) {
      message("The file '",paste0(basename(rm.folder), "_", names(res)[2],".csv"),"' has been generated and stored in the current working directory at '", getwd(),"'.")
      readr::write_delim(res, path = paste0(basename(rm.folder), "_", names(res)[2],".csv"), delim = ";")
    }

    return(res)
  }
}



