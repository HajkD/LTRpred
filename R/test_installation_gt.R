#' @title Helper function to test GenomeTools installation
#' @noRd
test_installation_gt <- function() {
    # test if a valid genometools (gt) version is installed
    tryCatch({
        sys_out <-
            system("gt suffixerator --version", intern = TRUE)
    }, error = function(e)
        stop(
            "It seems like you don't have genometools installed locally on your machine
            or the PATH variable to the genometools program is not set correctly.",
            "Please consult the Installation vignette or https://github.com/genometools/genometools 
            for details on how to install genometools.",
            call. = FALSE
        ))
    
    if (any(stringr::str_detect(sys_out, "GenomeTools")))
        return(TRUE)
    
}
