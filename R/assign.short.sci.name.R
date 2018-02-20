#' @title Transform long file names to short scientific names
#' @description This function takes a \code{tibble} as input
#' and transform long file names to short scientific names.
#' @param data a \code{tibble} from \code{\link{meta.seq.space}}.
#' @author Hajk-Georg Drost
assign.short.sci.name <- function(data) {
  new_spec_names <- unlist(lapply(data$organism, function(x) {
    split_str <- unlist(stringr::str_split(x, "_"))
    return(paste0(
      stringr::str_to_upper(stringr::str_sub(split_str[1], 1, 1)),
      split_str[2],
      collapse = ""
    ))
    
  }))
  
  data <- dplyr::mutate(data, organism = new_spec_names)
  return(data)
}
