#' @title Helper function to test Perl language installation
#' @noRd
test_installation_perl <- function() {
    # test if Perl language is installed
    tryCatch({
        sys_out <-
            system("perl -v", intern = TRUE)
    }, error = function(e)
        stop(
            "It seems like you don't have Perl installed locally on your machine
            or the PATH variable to Perl is not set correctly.",
            "Please consult the Installation vignette or https://www.perl.org/ 
            for details on how to install Perl on your system.",
            call. = FALSE
        ))
    
    if (any(stringr::str_detect(sys_out, "perl")))
        return(TRUE)
    
}
