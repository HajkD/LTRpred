#' @title Query the Dfam database to annotate putative LTRs
#' @description Validate or annotate putative LTR
#' transposons that have been predicted using \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}.
#' @param seq.file file path to the putative LTR transposon sequences in \code{fasta} format.
#' @param Dfam.db folder path to the local Dfam database.
#' @param eval E-value threshhold to perform the HMMer search against the Dfam database.
#' @param cores number of cores to use to perform parallel computations.
#' @param output.folder folder path to store the annotation results.
#' @author Hajk-Georg Drost
#' @details 
#' The Dfam database provides a collection of curated transposable element annotations.
#' @export
dfam.query <- function(seq.file, 
                       Dfam.db       = NULL, 
                       eval          = 1E-5, 
                       cores         = 1, 
                       output.folder = getwd()){
  
  if (is.null(Dfam.db))
    stop ("Please provide either a path to the Dfam.hmm database (file) or choose
          Dfam.db = 'download' so that Dfam.hmm is automatically loaded by this function
          (make sure that the internet connection is stabe.")
  
  if (Dfam.db == "download"){
    
    downloader::download("http://dfam.org/web_download/Current_Release/Dfam.hmm.gz",file.path(output.folder,"Dfam.hmm.gz"),mode = "wb")
    cat("\n")
    cat("Prepare the Dfam.hmm database...")
    cat("\n")
    system(paste0("hmmpress ",file.path(ws.wrap.path(output.folder),"Dfam.hmm")))
    cat("Run Dfam scan...")
    cat("\n")
    system(paste0("perl ",ws.wrap.path(system.file("dfamscan.pl", package = "LTRpred", mustWork = TRUE))," -fastafile ",ws.wrap.path(seq.file)," -hmmfile ",file.path(ws.wrap.path(output.folder),"Dfam.hmm")," -dfam_outfile ",
                  ws.wrap.path(file.path(output.folder,paste0(basename(seq.file),"_DfamAnnotation.out"))) ," -E ", eval,
                 " -cpu ",cores," --log_file ",ws.wrap.path(file.path(output.folder,paste0(basename(seq.file),"_logfile.txt"))) ," --masking_thresh "))
    cat("Finished Dfam scan!")
    cat("\n")
  } else {
    
    cat("\n")
    cat("Prepare the Dfam.hmm database...")
    cat("\n")
    system(paste0("hmmpress ",file.path(Dfam.db,"Dfam.hmm")))
    cat("Run Dfam scan...")
    cat("\n")
    system(paste0("perl ",ws.wrap.path(system.file("dfamscan.pl", package = "LTRpred", mustWork = TRUE))," -fastafile ",ws.wrap.path(seq.file)," -hmmfile ",ws.wrap.path(file.path(Dfam.db,"Dfam.hmm"))," -dfam_outfile ",
                  ws.wrap.path(file.path(output.folder,paste0(basename(seq.file),"_DfamAnnotation.out"))) ," -E ", eval,
                  " -cpu ",cores," --log_file ",ws.wrap.path(file.path(output.folder,paste0(basename(seq.file),"_logfile.txt"))) ," --masking_thresh "))
    cat("Finished Dfam scan!")
    cat("\n")
    cat("\n")
    }
  
}