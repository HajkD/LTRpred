#' @title Run LTRharvest to predict putative LTR Retrotransposons
#' @description This function implements an interface between R and
#' the LTRharvest command line tool to predict putative LTR retrotransposons from R.
#' @param genome.file path to the genome file in \code{fasta} format.
#' @param index.file specify the name of the enhanced suffix array index file that is computed
#'  by \code{suffixerator}. This opten can be used in case the suffix file was previously 
#'  generated, e.g. during a previous call of this function. In this case the suffix array index
#'  file does not need to be re-computed for new analyses. This is particularly useful when 
#'  running \code{LTRharvest} with different parameter settings.
#' @param range define the genomic interval in which predicted LTR transposons shall be reported
#' . In case \code{range[1] = 1000} and \code{range[2] = 10000} then candidates are only 
#' reported if they start after position 1000 and end before position 10000 in their respective 
#' sequence coordinates. If \code{range[1] = 0} and \code{range[2] = 0}, 
#' so \code{range = c(0,0)} (default) then the entire genome is being scanned.
#' @param seed  the minimum length for the exact maximal repeats. Only repeats with the specified minimum length are considered in all subsequent analyses. Default is \code{seed = 30}.
#' @param minlenltr minimum LTR length. Default is \code{minlenltr = 100}. 
#' @param maxlenltr maximum LTR length. Default is \code{maxlenltr = 3500}.
#' @param mindistltr minimum distance of LTR starting positions. Default is \code{mindistltr = 4000}.
#' @param maxdistltr maximum distance of LTR starting positions. Default is \code{maxdistltr = 25000}.
#' @param similar minimum similarity value between the two LTRs in percent. \code{similar = 70}.
#' @param mintsd minimum target site duplications (TSDs) length. If no search for TSDs
#' shall be performed, then specify \code{mintsd = NULL}. Default is \code{mintsd = 4}.
#' @param maxtsd maximum target site duplications (TSDs) length. If no search for TSDs
#' shall be performed, then specify \code{maxtsd = NULL}. Default is \code{maxtsd = 20}.
#' @param vic number of nucleotide positions left and right (the vicinity) of the predicted
#'  boundary of a LTR that will be searched for TSDs and/or one motif (if specified). 
#'  Default is \code{vic = 60}.
#' @param overlaps specify how overlapping LTR retrotransposon predictions shall be treated. 
#' If \code{overlaps = "no"} is selected, then neither nested nor overlapping predictions will be reported in the output. In case \code{overlaps = "best"} is selected then in the case of two or more nested or overlapping predictions, solely the LTR retrotransposon prediction with
#' the highest similarity between its LTRs will be reported.
#' If \code{overlaps = "all"} is selected then all LTR retrotransposon predictions 
#' will be reported whether there are nested and/or overlapping predictions or not. 
#' Default is \code{overlaps = "best"}.
#' @param xdrop specify the xdrop value (> 0) for extending a seed repeat in both directions
#'  allowing for matches, mismatches, insertions, and deletions. The xdrop extension process
#'   stops as soon as the extension involving matches, mismatches, insersions, and deletions 
#'   has a score smaller than T -X, where T denotes the largest score seen so far. Default is \code{cdrop = 5}.
#' @param mat specify the positive match score for the X-drop extension process. Default is \code{mat = 2}.
#' @param mis specify the negative mismatch score for the X-drop extension process. Default is \code{mis = -2}.
#' @param ins specify the negative insertion score for the X-drop extension process. Default is \code{ins = -3}.
#' @param del specify the negative deletion score for the X-drop extension process. Default is \code{del = -3}.
#' @param motif specify 2 nucleotides for the starting motif and 2 nucleotides for the ending
#'  motif at the beginning and the ending of each LTR, respectively.
#'  Only palindromic motif sequences - where the motif sequence is equal to its complementary
#'  sequence read backwards - are allowed, e.g. \code{motif = "tgca"}. Type the nucleotides without any space
#'  separating them. If this option is not selected by the user, candidate pairs will not be
#'  screened for potential motifs. If this options is set but no allowed number of
#'  mismatches is specified by the argument \code{motifmis} and a search for the exact 
#'  motif will be conducted. If \code{motif = NULL} then no explicit motif is being specified.
#' @param motifmis allowed number of mismatches in the TSD motif specified in \code{motif}. 
#' The number of mismatches needs to be between [0,3].  Default is \code{motifmis = 0}.
#' @param output.path a path/folder to store all results returned by \code{LTRharvest}. 
#' If \code{output.path = NULL} (Default) then a folder with the name of the input genome file
#' will be generated in the current working directory of R and all results are then stored in this folder.
#' @param verbose logical value indicating whether or not detailed information shall be printed on the console.
#' @author Hajk-Georg Drost
#' @details 
#' The \code{LTRharvest} function provides an interface to the \code{LTRharvest} command line
#' tool and furthermore takes care of the entire folder handling, output parsing, and data
#' processing of the \code{LTRharvest} prediction.
#' 
#' Internally a folder named \code{output.path}_ltrharvest is generated and all computations
#' returned by \code{LTRharvest} are then stored in this folder. These files (see section \code{Value}) are then parsed and returned as list of data.frames by this function.
#' 
#' \code{LTRharvest} can be used as independently or as initial pre-computation step
#' to sufficiently detect LTR retrotransposons with \code{LTRdigest}. 
#' @examples 
#' \dontrun{
#' 
#' # Run LTRharvest for Arabidopsis thaliana using standard parameters
#' LTRharvest(genome.file = "TAIR10_chr_all.fas")
#' 
#' # Run LTRharvest for Arabidopsis thaliana using standard parameters
#' # and use an existing (already computed) suffixarray index file of 
#' # the A. thaliana genome
#' LTRharvest(genome.file      = "TAIR10_chr_all.fas", 
#'            index.file       = "TAIR10_chr_all_index.fas")
#' }
#' @references 
#' D Ellinghaus, S Kurtz and U Willhoeft. LTRharvest, an efficient and flexible software for de novo detection of LTR retrotransposons. BMC Bioinformatics (2008). 9:18.
#' 
#' Most argument specifications are adapted from the User manual of LTRharvest.
#' @seealso \code{\link{LTRdigest}},  \code{\link{LTRpred}}, \code{\link{PlotLTRAge}}, 
#' \code{\link{PlotLTRWidth}}, \code{\link{PlotLTRRange}},
#' \code{\link{read.prediction}}, \code{\link{read.seqs}},
#' \code{\link{pred2fasta}}, \code{\link{pred2gff}}
#' @return
#' The \code{LTRharvest} function generates the following output files:
#' 
#' \itemize{
#' \item *_BetweenLTRSeqs.fsa : DNA sequences of the region between the LTRs in fasta format. 
#' \item *_Details.tsv : A spread sheet containing detailed information about the predicted LTRs.
#' \item *_FullLTRRetrotransposonSeqs.fsa : DNA sequences of the entire predicted LTR retrotransposon.
#' \item *_index.fsa : The suffixarray index file used to predict putative LTR retrotransposonswith \code{LTRharvest}.
#' \item *_Prediction.gff : A spread sheet containing detailed additional information about the predicted LTRs (partially redundant with the *_Details.tsv file).
#' }
#' The ' * ' is an place holder for the name of the input genome file.
#' @export

LTRharvest <- function(genome.file,
                       index.file  = NULL,
                       range       = c(0,0),
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
                       motif       = NULL,
                       motifmis    = 0,
                       output.path = NULL,
                       verbose     = TRUE){
    
    
  if ((is.null(mintsd) & !is.null(maxtsd)) || (!is.null(mintsd) & is.null(maxtsd)))
    stop ("Please specify a TSD length for both: min and max TSD length!")
  
  if (!is.element(overlaps,c("no","best","all")))
    stop ("Please select as LTR transposon overlap option either, 'no', 'best', or 'all'.")
    
    if (is.null(output.path)){
        output.path <- paste0(unlist(stringr::str_split(basename(genome.file),"[.]"))[1],"_ltrharvest")
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
    IndexOutputFileName <- file.path(output.path,paste0(OutputFileNameIdentifier,"_index",".fsa"))
    
    cat("\n")
    cat("Run LTRharvest...")
    cat("\n")
    
    # In case no index file is passed to this function
    # an index file will be generated using "gt suffixerator"
    if (is.null(index.file)){
        if (verbose){
          cat("\n")
          cat("Generating the indexfile ",IndexOutputFileName," with suffixerator...")
          cat("\n")
        }
        # Genrate Suffix for LTRharvest
        system(paste0("gt suffixerator -db ",ws.wrap.path(genome.file)," -indexname ", IndexOutputFileName," -tis -suf -lcp -des -ssp -sds -dna"))    
    } else {
        IndexOutputFileName <- index.file
    }
    
    if (verbose){
      cat("Running LTRharvest and write results to ",output.path,"...")
      cat("\n")  
    }
      
    if (is.null(motif)){
        
        # Run LTRharvest without motif
        system(paste0("gt ltrharvest -index ", ws.wrap.path(IndexOutputFileName)," \ ",
                      "-range ",range[1]," ",range[2]," \ ",
                      "-seed ", seed," \ ",
                      "-minlenltr ", minlenltr, " \ ",
                      "-maxlenltr ", maxlenltr, " \ ",
                      "-mindistltr ", mindistltr, " \ ",
                      "-maxdistltr ", maxdistltr, " \ ",
                      "-similar ", similar, " \ ",
                      ifelse(!is.null(mintsd),paste0("-mintsd ", mintsd, " \ ")," "),
                      ifelse(!is.null(maxtsd),paste0("-maxtsd ", maxtsd, " \ ")," "),
                      "-vic ", vic, " \ ",
                      "-overlaps ", overlaps, " \ ",
                      "-xdrop ", xdrop, " \ ",
                      "-mat ", mat, " \ ",
                      "-mis ", mis, " \ ",
                      "-ins ", ins, " \ ",
                      "-del ", del, " \ ",
                      ifelse(!is.null(mintsd),"-longoutput \ "," "),
                      "-out ", ws.wrap.path(file.path(output.path,paste0(OutputFileNameIdentifier,"_FullLTRretrotransposonSeqs",".fsa"))) , " \ ",
                      "-outinner ", ws.wrap.path(file.path(output.path,paste0(OutputFileNameIdentifier,"_BetweenLTRSeqs",".fsa"))) , " \ ",
                      "-gff3 ", ws.wrap.path(file.path(output.path,paste0(OutputFileNameIdentifier,"_Prediction",".gff"))), " \ ",
                      ">> ", file.path(output.path,paste0(OutputFileNameIdentifier,"_Details",".tsv"))," 2>&1"))
    }
    
    if (!is.null(motif)){
        
        if (!is.element(motifmis,seq_len(3)))
            stop ("The 'motifmis' argument should be between [0,3].")
        
        if (nchar(motif) > 4)
            stop ("Please choose 2 nucleotides for the starting motif and 2 nucleotides for the ending motif.")
        
        # Run LTRharvest with motif
        system(paste0("gt ltrharvest -index ", ws.wrap.path(IndexOutputFileName)," \ ",
                      "-seed ", seed," \ ",
                      "-minlenltr ", minlenltr, " \ ",
                      "-maxlenltr ", maxlenltr, " \ ",
                      "-mindistltr ", mindistltr, " \ ",
                      "-maxdistltr ", maxdistltr, " \ ",
                      "-similar ", similar, " \ ",
                      ifelse(!is.null(mintsd),paste0("-mintsd ", mintsd, " \ ")," "),
                      ifelse(!is.null(maxtsd),paste0("-maxtsd ", maxtsd, " \ ")," "),
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
                      "-out ", ws.wrap.path(file.path(output.path,paste0(OutputFileNameIdentifier,"_FullLTRretrotransposonSeqs",".fsa"))) , " \ ",
                      "-outinner ", ws.wrap.path(file.path(output.path,paste0(OutputFileNameIdentifier,"_BetweenLTRSeqs",".fsa"))) , " \ ",
                      "-gff3 ", ws.wrap.path(file.path(output.path,paste0(OutputFileNameIdentifier,"_Prediction",".gff")))))
        
    }
    
    if (verbose){
      cat("LTRharvest analysis finished!")
      cat("\n")  
    }
}

