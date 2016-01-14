#' @title Run LTRdigest to predict putative LTR Retrotransposons
#' @description This function implements an interface between R and
#' the LTRdigest command line tool to predict putative LTR retrotransposons from R.
#' @param input.gff3
#' @param genome.file
#' @param aaout
#' @param aliout
#' @param pptlen
#' @param uboxlen
#' @param pptradius
#' @param trnas
#' @param pbsalilen
#' @param pbsoffset
#' @param pbstrnaoffset
#' @param pbsmaxedist
#' @param pbsradius
#' @param hmms
#' @param pdomevalcutoff
#' @param pbsmatchscore
#' @param pbsmismatchscore
#' @param pbsinsertionscore
#' @param pbsdeletionscore
#' @param pfam.ids
#' @param cores
#' @param index.file
#' @param output.path
#' @author Hajk-Georg Drost
#' @examples 
#' \dontrun{
#' 
#' }
#' @export

LTRdigest <- function(input.gff3,
                      genome.file,
                      aaout             = "yes",
                      aliout            = "yes",
                      pptlen            = c(8,30),
                      uboxlen           = c(3,30),
                      pptradius         = 30,
                      trnas             = NULL,
                      pbsalilen         = c(11,30),
                      pbsoffset         = c(0,5),
                      pbstrnaoffset     = c(0,5),
                      pbsmaxedist       = 1,
                      pbsradius         = 30,
                      hmms              = NULL,
                      pdomevalcutoff    = 1E-5,
                      pbsmatchscore     = 5,
                      pbsmismatchscore  = -10,
                      pbsinsertionscore = -20,
                      pbsdeletionscore  = -20,
                      pfam.ids          = NULL,
                      cores             = 1,
                      index.file        = NULL,
                      output.path       = NULL){
    
    
    if (parallel::detectCores() < cores)
        stop ("Your system does not provide the number of cores you specified.")
    
    if (is.null(output.path)){
        output.path <- paste0(unlist(stringr::str_split(basename(genome.file),"[.]"))[1],"_ltrdigest")
        if (dir.exists(output.path)){
            unlink(output.path, recursive = TRUE)
            dir.create(output.path)
        } else {
            dir.create(output.path)
        }
    } else {
        
        if (dir.exists(output.path)){
            unlink(output.path, recursive = TRUE)
            dir.create(output.path)
        } else {
            dir.create(output.path)
        }
    }
    
    OutputFileNameIdentifier <- unlist(stringr::str_split(basename(genome.file),"[.]"))[1] 
    IndexOutputFileName <- file.path(output.path,paste0(OutputFileNameIdentifier,"_index_ltrdigest",".fsa"))
    
    if (!is.null(pfam.ids)){
        
        if (!is.null(hmms))
            stop ("Please only provide the PFAM ids and not the actual hmm files when using the 'pfam.ids' argument.")
        
        sapply(pfam.ids, function(pfamid) utils::download.file(url      = paste0("http://pfam.xfam.org/family/",pfamid,"/hmm"),
                                                               destfile = paste0("hmm_",pfamid,".hmm"),
                                                               quiet    = TRUE))
        
    }
    
    
    # In case no index file is passed to this function
    # an index file will be generated using "gt suffixerator"
    if (is.null(index.file)){
        cat("\n")
        cat("Generating the index file ",IndexOutputFileName," with suffixerator...")
        cat("\n")
        # Genrate Suffix for LTRdigest
        system(paste0("gt suffixerator -tis -des -dna -ssp -db ",genome.file," -indexname ", IndexOutputFileName))    
    } else {
        IndexOutputFileName <- index.file
    }
    
    sorted.input.gff3 <- paste0(unlist(stringr::str_split(input.gff3,"[.]"))[1],"_sorted.gff")
    
    cat("Sort index file...")
    cat("\n")
    system(paste0("gt gff3 -sort ",input.gff3," > ", sorted.input.gff3))
    
    
    cat("Running LTRdigest and write results to ",output.path," ...")
    cat("\n")    
    
    
        # Run LTRdigest
        system(paste0("gt -j ",cores," ltrdigest "," \ ",
                      "-aaout ", aaout," \ ",
                      "-aliout ", aliout," \ ",
                      "-pptlen ", pptlen[1]," ",pptlen[2]," \ ",
                      "-uboxlen ", uboxlen[1]," ",uboxlen[2], " \ ",
                      "-pptradius ", pptradius, " \ ",
                      ifelse(is.null(trnas),"", paste0("-trnas ", trnas, " \ ")),
                      "-pbsalilen ", pbsalilen[1]," ",pbsalilen[2], " \ ",
                      "-pbsoffset ", pbsoffset[1]," ",pbsoffset[2], " \ ",
                      "-pbstrnaoffset ", pbstrnaoffset[1]," ",pbstrnaoffset[2], " \ ",
                      "-pbsmaxedist ", pbsmaxedist, " \ ",
                      "-pbsradius ", pbsradius, " \ ",
                      ifelse(is.null(hmms),ifelse(!is.null(pfam.ids),"-hmms hmm_*",""), paste0("-hmms ", paste0(hmms,collapse = " "), " \ ")),
                      "-pdomevalcutoff ", pdomevalcutoff, " \ ",
                      "-pbsmatchscore ", pbsmatchscore, " \ ",
                      "-pbsmismatchscore ", pbsmismatchscore, " \ ",
                      "-pbsinsertionscore ", pbsinsertionscore, " \ ",
                      "-pbsdeletionscore ", pbsdeletionscore, " \ ",
                      "-outfileprefix ", file.path(output.path,paste0(OutputFileNameIdentifier,"-ltrdigest")) ," ",sorted.input.gff3," \ "
                      ,IndexOutputFileName, " > ", file.path(output.path,paste0(OutputFileNameIdentifier,"_LTRdigestPrediction",".gff"))))
        
    
    cat("Analysis finished!")
    cat("\n") 
    
    
}
