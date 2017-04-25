#' @title Helper function to test HMMer installation
#' @noRd
test_installation_hmmer <- function() {
    # test if a valid HMMer version is installed
    tryCatch({
        sys_out <-
            system("hmmpress -h", intern = TRUE)
    }, error = function(e)
        stop(
            "It seems like you don't have HMMer installed locally on your machine
            or the PATH variable to the HMMer program is not set correctly.",
            "Please consult the Installation vignette or http://hmmer.org/download.html
            for details on how to install HMMer",
            call. = FALSE
        ))
    
    if (any(stringr::str_detect(sys_out, "HMMER")))
        return(TRUE)
    
}
