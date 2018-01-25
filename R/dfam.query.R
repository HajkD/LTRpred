#' @title Query the Dfam database to annotate putative LTRs
#' @description Validate or annotate putative LTR
#' transposons that have been predicted using \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}.
#' @param seq.file file path to the putative LTR transposon sequences in \code{fasta} format.
#' @param Dfam.db folder path to the local Dfam database or \code{Dfam.db = "download"} in case the Dfam
#'  database shall be automatically downloaded before performing query analyses.
#' @param eval E-value threshhold to perform the HMMer search against the Dfam database.
#' @param cores number of cores to use to perform parallel computations.
#' @param output.folder folder path to store the annotation results.
#' @author Hajk-Georg Drost
#' @details 
#' The Dfam database provides a collection of curated transposable element annotations.
#' @export
dfam.query <- function(seq.file, 
                       Dfam.db       = NULL, 
                       eval          = 1E-3, 
                       cores         = 1, 
                       output.folder = getwd()){
  
    
    test_installation_perl()
    test_installation_hmmer()
    
    if (is.null(Dfam.db))
        stop(
            "Please provide either a path to the Dfam.hmm database (file) or choose
            Dfam.db = 'download' so that Dfam.hmm is automatically loaded by this function
            (make sure that the internet connection is stabe."
             )
    
    if (!file.exists("/usr/local/bin/dfamscan.pl"))
        stop(
            "The perl script 'dfamscan.pl' could not be found! Please download 'dfamscan.pl' from www.dfam.org/web_download/Current_Release/dfamscan.pl and store it in '/usr/local/bin'."
        )
    
    if (!file.exists(seq.file))
        stop("The file '",seq.file,"' does not exist! Please make sure that a correct path to the sequence file is passed to the dfam.query() function.", call. = FALSE)
    
    if (Dfam.db == "download") {
      message("Download Dfam database from http://dfam.org/web_download/Current_Release/Dfam.hmm.gz ...")
        downloader::download(
            "http://dfam.org/web_download/Current_Release/Dfam.hmm.gz",
            file.path(output.folder, "Dfam.hmm.gz"),
            mode = "wb"
        )
        message("Download completed!")
        message("Prepare the Dfam.hmm database...")
        
        hmmpress <- system(paste0("hmmpress ", file.path(
            ws.wrap.path(output.folder), gzfile("Dfam.hmm.gz")
        )), intern = TRUE)
        
        if (attr(hmmpress, "status") > 0)
          stop("hmmpress could not format the file ",file.path(
            ws.wrap.path(output.folder), gzfile("Dfam.hmm.gz")), ". Is hmmpress installed on your system and did the download process of the Dfam database work properly? ", call. = FALSE)
        
        message("Run Dfam scan...")
        # make sure that in future versions the PATH variable is set and OSX, Linux,
        # and Windows paths will be supported
        dfamscan <- system(
            paste0(
                "perl ",
                "/usr/local/bin/dfamscan.pl",
                " -fastafile ",
                ws.wrap.path(seq.file),
                " -hmmfile ",
                file.path(ws.wrap.path(output.folder), "Dfam.hmm"),
                " -dfam_outfile ",
                ws.wrap.path(file.path(
                    output.folder, paste0(basename(seq.file), "_DfamAnnotation.out")
                )) ,
                " -E ",
                eval,
                " -cpu ",
                cores,
                " --log_file ",
                ws.wrap.path(file.path(
                    output.folder, paste0(basename(seq.file), "_logfile.txt")
                )) ,
                " --masking_thresh "
            ), intern = TRUE
        )
        
        if (attr(dfamscan, "status") > 0)
          stop("dfamscan could not be run properly! Did you store the file 'dfamscan.pl' in '/usr/local/bin/dfamscan.pl' ?", call. = FALSE)
        
        message("Finished Dfam scan!")
        message("A dfam query file has been generated and stored at",ws.wrap.path(file.path(
            output.folder, paste0(basename(seq.file), "_DfamAnnotation.out"))),".")
    } else {
        if (file.exists(file.path(ws.wrap.path(output.folder), "Dfam.hmm.h3f")) &
            file.exists(file.path(ws.wrap.path(output.folder), "Dfam.hmm.h3i")) &
            file.exists(file.path(ws.wrap.path(output.folder), "Dfam.hmm.h3m")) &
            file.exists(file.path(ws.wrap.path(output.folder), "Dfam.hmm.h3p"))) {
            message("Prepare the Dfam.hmm database...")
            system(paste0("hmmpress ", file.path(Dfam.db, "Dfam.hmm")))
        }
        
        message("Run Dfam scan...")
        system(
            paste0(
                "perl ",
                "/usr/local/bin/dfamscan.pl",
                " -fastafile ",
                ws.wrap.path(seq.file),
                " -hmmfile ",
                ws.wrap.path(file.path(Dfam.db, "Dfam.hmm")),
                " -dfam_outfile ",
                ws.wrap.path(file.path(
                    output.folder, paste0(basename(seq.file), "_DfamAnnotation.out")
                )) ,
                " -E ",
                eval,
                " -cpu ",
                cores,
                " --log_file ",
                ws.wrap.path(file.path(
                    output.folder, paste0(basename(seq.file), "_logfile.txt")
                )) ,
                " --masking_thresh "
            )
        )
        message("Finished Dfam scan!")
        message("A dfam query file has been generated and stored at",ws.wrap.path(file.path(
            output.folder, paste0(basename(seq.file), "_DfamAnnotation.out"))),".")
        
    }
}
