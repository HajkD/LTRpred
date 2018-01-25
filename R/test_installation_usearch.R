#' @title Helper function to test USEARCH installation
#' @noRd
test_installation_usearch <- function() {
    # test if a valid USEARCH version is installed
    tryCatch({
        sys_out <-
            system("usearch --version", intern = TRUE)
    }, error = function(e)
        stop(
            "It seems like you don't have USEARCH installed locally on your machine or the PATH variable to the USEARCH program is not set correctly. ",
            "Please consult the Installation vignette or http://drive5.com/usearch/download.html for details on how to install USEARCH.",
            call. = FALSE
        ))
    
    if (any(stringr::str_detect(sys_out, "usearch")))
        return(TRUE)
    
}
