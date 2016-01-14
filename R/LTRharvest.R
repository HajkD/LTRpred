#' @title Run LTRharvest to predict putative LTR Retrotransposons
#' @description This function implements an interface between R and
#' the LTRharvest command line tool to predict putative LTR retrotransposons from R.
#' @param input
#' @param index.file
#' @param seed
#' @param minlenltr
#' @param maxlenltr
#' @param mindistltr
#' @param maxdistltr
#' @param similar
#' @param mintsd
#' @param maxtsd
#' @param vic
#' @param overlaps
#' @param xdrop
#' @param mat
#' @param mis
#' @param ins
#' @param del
#' @param motifmis
#' @param motif
#' @param output.path
#' @author Hajk-Georg Drost
#' @details 
#' 
#' @examples 
#' \dontrun{
#' 
#' }
#' 
#' @export

LTRharvest <- function(input,
                       index.file  = NULL,
                       seed        = 30,
                       minlenltr   = 100,
                       maxlenltr   = 3500,
                       mindistltr  = 4000,
                       maxdistltr  = 25000,
                       similar     = 70,
                       mintsd      = 4,
                       maxtsd      = 20,
                       vic         = 60,
                       overlaps    = "no",
                       xdrop       = 5,
                       mat         = 2,
                       mis         = -2,
                       ins         = -3,
                       del         = -3,
                       motifmis    = 0,
                       motif       = NULL,
                       output.path = NULL){
    
    
    if (is.null(output.path)){
        output.path <- unlist(stringr::str_split(basename(input),"[.]"))[1]
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

    OutputFileNameIdentifier <- unlist(stringr::str_split(basename(input),"[.]"))[1] 
    IndexOutputFileName <- file.path(output.path,paste0(OutputFileNameIdentifier,"_index",".fsa"))
    
    # In case no index file is passed to this function
    # an index file will be generated using "gt suffixerator"
    if (is.null(index.file)){
        cat("Generating the indexfile ",IndexOutputFileName," with suffixerator...")
        cat("\n")
        # Genrate Suffix for LTRharvest
        system(paste0("gt suffixerator -db ",input," -indexname ", IndexOutputFileName," -tis -suf -lcp -des -ssp -sds -dna"))    
    } else {
        IndexOutputFileName <- index.file
    }
    
    cat("\n")
    cat("Running LTRharvest and write results to ",output.path," ...")
    cat("\n")    
    
    if (is.null(motif)){
        
        # Run LTRharvest without motif
        system(paste0("gt ltrharvest -index ", IndexOutputFileName," \ ",
                      "-seed ", seed," \ ",
                      "-minlenltr ", minlenltr, " \ ",
                      "-maxlenltr ", maxlenltr, " \ ",
                      "-mindistltr ", mindistltr, " \ ",
                      "-maxdistltr ", maxdistltr, " \ ",
                      "-similar ", similar, " \ ",
                      "-mintsd ", mintsd, " \ ",
                      "-maxtsd ", maxtsd, " \ ",
                      "-vic ", vic, " \ ",
                      "-overlaps ", overlaps, " \ ",
                      "-xdrop ", xdrop, " \ ",
                      "-mat ", mat, " \ ",
                      "-mis ", mis, " \ ",
                      "-ins ", ins, " \ ",
                      "-del ", del, " \ ",
                      "-longoutput ", " \ ",
                      "-out ", file.path(output.path,paste0(OutputFileNameIdentifier,"_FullLTRretrotransposonSeqs",".fsa")) , " \ ",
                      "-outinner ", file.path(output.path,paste0(OutputFileNameIdentifier,"_BetweenLTRSeqs",".fsa")) , " \ ",
                      "-gff3 ", file.path(output.path,paste0(OutputFileNameIdentifier,"_Prediction",".gff")), " \ ",
                      ">> ", file.path(output.path,paste0(OutputFileNameIdentifier,"_Details",".tsv"))," 2>&1"))
    }
    
    if (!is.null(motif)){
        
        if (!is.element(motifmis,seq_len(3)))
            stop ("The 'motifmis' argument should be between [0,3].")
        
        if (nchar(motif) > 4)
            stop ("Please choose 2 nucleotides for the starting motif and 2 nucleotides for the ending motif.")
        
        # Run LTRharvest with motif
        system(paste0("gt ltrharvest -index ", IndexOutputFileName," \ ",
                      "-seed ", seed," \ ",
                      "-minlenltr ", minlenltr, " \ ",
                      "-maxlenltr ", maxlenltr, " \ ",
                      "-mindistltr ", mindistltr, " \ ",
                      "-maxdistltr ", maxdistltr, " \ ",
                      "-similar ", similar, " \ ",
                      "-mintsd ", mintsd, " \ ",
                      "-maxtsd ", maxtsd, " \ ",
                      "-motif ", motif, " \ ",
                      "-motifmis ", motifmis, " \ ",
                      "-vic ", vic, " \ ",
                      "-overlaps ", overlaps, " \ ",
                      "-xdrop ", xdrop, " \ ",
                      "-mat ", mat, " \ ",
                      "-mis ", mis, " \ ",
                      "-ins ", ins, " \ ",
                      "-del ", del, " \ ",
                      "-longoutput ", " \ ",
                      "-out ", file.path(output.path,paste0(OutputFileNameIdentifier,"_FullLTRretrotransposonSeqs",".fsa")) , " \ ",
                      "-outinner ", file.path(output.path,paste0(OutputFileNameIdentifier,"_BetweenLTRSeqs",".fsa")) , " \ ",
                      "-gff3 ", file.path(output.path,paste0(OutputFileNameIdentifier,"_Prediction",".gff"))))
        
    }
    
    cat("Analysis finished!")
    cat("\n") 
    
    
}

