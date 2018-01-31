str_chop_vec <- function(x, pattern) {
    return(unlist(lapply(x, function(y) {
        
        unlist(stringr::str_split(y, pattern))[1]
        
    })))
}

