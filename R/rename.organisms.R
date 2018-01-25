#' @title Rename file path names from \code{read_proteome}, etc.
#' @description Rename a file path e.g. \code{Arabidopsis_thaliana.fa.gz} to \code{Athaliana}.
#' @param data a GenomeInfo file.
#' @author Hajk-Georg Drost
#' @export 
rename.organisms <- function(data) {
    new_spec_names <- unlist(lapply(data$organism, function(x) {
        split_str <- unlist(stringr::str_split(x, "_"))
        return(paste0(
            stringr::str_to_upper(stringr::str_sub(split_str[1], 1, 1)),
            split_str[2],
            collapse = ""
        ))
        
    }))
    
    data <- dplyr::mutate(data, organism = new_spec_names)
    
}
