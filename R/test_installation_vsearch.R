#' @title Helper function to test VSEARCH installation
#' @noRd
test_installation_vsearch <- function() {
    # test if a valid VSEARCH version is installed
    tryCatch({
        sys_out <-
            system("vsearch -v", intern = TRUE)
    }, error = function(e)
        stop(
            "It seems like you don't have VSEARCH installed locally on your machine
            or the PATH variable to the VSEARCH program is not set correctly.",
            "Please consult the Installation vignette or https://github.com/torognes/vsearch 
            for details on how to install VSEARCH.",
            call. = FALSE
        ))
    
    if (any(stringr::str_detect(sys_out, "vsearch")))
        return(TRUE)
    
}
