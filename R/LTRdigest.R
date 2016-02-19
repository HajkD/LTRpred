#' @title Run LTRdigest to predict putative LTR Retrotransposons
#' @description This function implements an interface between R and
#' the LTRdigest command line tool to predict putative LTR retrotransposons from R.
#' @param input.gff3 path to the prediction file in gff3 format returned by \code{\link{LTRharvest}}.
#' @param genome.file path to the genome file in \code{fasta} format.
#' @param aaout shall the protein sequence of the HMM matches to the predicted LTR transposon 
#' be generated as fasta file or not. Options are \code{aaout = "yes"} or \code{aaout = "no"}.
#' @param aliout shall the alignment of the protein sequence of the HMM matches to the predicted LTR transposon 
#' be generated as fasta file or not. Options are \code{aaout = "yes"} or \code{aaout = "no"}.
#' @param pptlen a two dimensional numeric vector specifying the minimum and maximum allowed
#' lengths for PPT predictions. If a purine-rich region that does not fulfill this range is
#' found, it will be discarded. Default is \code{pptlen = c(8,30)} (minimum = 8; maximum = 30).
#' @param uboxlen a two dimensional numeric vector specifying the minimum and maximum allowed
#' lengths for U-box predictions. If a T-rich region preceding a PPT that does not fulfill the PPT length criteria is
#' found, it will be discarded. Default is \code{uboxlen = c(3,30)} (minimum = 3; maximum = 30).
#' @param pptradius a numeric value specifying the area around the 3' LTR beginning to be 
#' considered when searching for PPT. Default value is \code{pptradius = 30}.
#' @param trnas path to the fasta file storing the unique tRNA sequences that shall be matched to the
#' predicted LTR transposon (tRNA library). 
#' @param pbsalilen a two dimensional numeric vector specifying the minimum and maximum allowed
#' lengths for PBS/tRNA alignments. If the local alignments are shorter or longer than this
#' range, it will be discarded. Default is \code{pbsalilen = c(11,30)} (minimum = 11; maximum = 30).
#' @param pbsoffset a two dimensional numeric vector specifying the minimum and maximum allowed
#' distance between the start of the PBS and the 3' end of the 5' LTR. Local alignments not 
#' fulfilling this criteria will be discarded. Default is \code{pbsoffset = c(0,5)} (minimum = 0; maximum = 5).
#' @param pbstrnaoffset a two dimensional numeric vector specifying the minimum and maximum allowed
#' PBS/tRNA alignment offset from the 3' end of the tRNA. Local alignments not 
#' fulfilling this criteria will be discarded. Default is \code{pbstrnaoffset = c(0,5)} (minimum = 0; maximum = 5).
#' @param pbsmaxedist a numeric value specifying the maximal allowed unit edit distance in a
#' local PBS/tRNA alignment.
#' @param pbsradius a numeric value specifying the area around the 5' LTR end to be 
#' considered when searching for PBS Default value is \code{pbsradius = 30}.
#' @param hmms a character string or a character vector storing either the hmm files for
#' searching internal domains between the LTRs of predicted LTR transposons or a vector of
#' Pfam IDs from http://pfam.xfam.org/ that are downloaded and used to search for corresponding protein domains
#' within the predicted LTR transposons. As an option users can rename all of their hmm files
#' so that they start for example with the name \code{hmms = "hmm_*"}. This way all files starting with 
#' \code{hmm_} will be considered for the subsequent protein domain search. In case Pfam IDs 
#' are specified, the \code{LTRpred} function will automatically download the corresponding 
#' HMM files and use them for further protein domain searches. In case users prefer to specify 
#' Pfam IDs please specify them in the \code{pfam.ids} parmeter and choose \code{hmms = NULL}.  
#' @param pdomevalcutoff a numeric value specifying the E-value cutoff for corresponding HMMER searches. All hits that do not fulfill this criteria are discarded. Default is \code{pdomevalcutoff = 1E-5}.
#' @param pbsmatchscore specify the match score used in the PBS/tRNA Smith-Waterman alignment.
#' Default is \code{pbsmatchscore = 5}.
#' @param pbsmismatchscore specify the mismatch score used in the PBS/tRNA Smith-Waterman alignment.
#' Default is \code{pbsmismatchscore = -10}.
#' @param pbsinsertionscore specify the insertion score used in the PBS/tRNA Smith-Waterman alignment.
#' Default is \code{pbsinsertionscore = -20}.
#' @param pbsdeletionscore specify the deletion score used in the PBS/tRNA Smith-Waterman alignment.
#' Default is \code{pbsdeletionscore = -20}.
#' @param pfam.ids a character vector storing the Pfam IDs from http://pfam.xfam.org/
#' that shall be downloaded and used to perform protein domain searches within the sequences
#' between the predicted LTRs.
#' @param cores number of cores to be used for multicore processing.
#' @param index.file specify the name of the enhanced suffix array index file that is computed
#'  by \code{suffixerator}. This opten can be used in case the suffix file was previously 
#'  generated, e.g. during a previous call of this function. In this case the suffix array index
#'  file does not need to be re-computed for new analyses. This is particularly useful when 
#'  running \code{LTRdigest} with different parameter settings.
#' @param output.path a path/folder to store all results returned by \code{LTRdigest}. 
#' If \code{output.path = NULL} (Default) then a folder with the name of the input genome file
#' will be generated in the current working directory of R and all results are then stored in this folder.
#' @author Hajk-Georg Drost
#' @details The \code{LTRdigest} function is a wrapper function to work with the
#' call the \code{LTRdigest} command line tool from R.
#' @examples 
#' \dontrun{
#' # Run LTRharvest for Arabidopsis thaliana using standard parameters
#' LTRharvest(genome.file = "TAIR10_chr_all.fas")
#' 
#' # Run LTRdigest for Arabidopsis thaliana using standard parameters
#' LTRdigest()
#' }
#' @return 
#' The \code{LTRdigest} function generates the following output files:
#' 
#' \itemize{
#' \item *_index_ltrdigest.fsa : The suffixarray index file used to predict putative LTR retrotransposonswith \code{LTRdigest}.
#' \item *_LTRdigestPrediction.gff : A spread sheet containing detailed information about the predicted LTRs.
#' \item *-ltrdigest_tabout.csv : A spread sheet containing additional detailed information about the predicted LTRs.
#' \item *-ltrdigest_complete.fas : The full DNA sequences of all predicted LTR transposons.
#' \item *-ltrdigest_conditions.csv : Contains information about the parameters used for a given
#' \code{LTRdigest} run.
#' \item *-ltrdigest_pbs.fas : Stores the predicted PBS sequences for the putative LTR retrotransposons.
#' \item *-ltrdigest_ppt.fas : Stores the predicted PPT sequences for the putative LTR retrotransposons.
#' \item *-ltrdigest_5ltr.fas and *-ltrdigest_3ltr.fas: Stores the predicted 5' and 3' LTR sequences. Note: If the direction of the putative retrotransposon could be predicted, these
#' files will contain the corresponding 3' and 5' LTR sequences. If no direction could be
#' predicted, forward direction with regard to the original sequence will be assumed by 
#' \code{LTRdigest}, i.e. the 'left' LTR will be considered the 5' LTR.
#' \item *-ltrdigest_pdom_<domainname>.fas : Stores the DNA sequences of the HMM matches
#' to the LTR retrotransposon candidates. 
#' \item *-ltrdigest_pdom_<domainname>_aa.fas : Stores the concatenated protein sequences of 
#' the HMM matches to the LTR retrotransposon candidates.
#' \item *-ltrdigest_pdom_<domainname>_ali.fas : Stores the alignment information for all matches of the given protein domain model to the translations of all candidates.
#' }
#' The ' * ' is an place holder for the name of the input genome file.
#' @seealso \code{\link{LTRharvest}},  \code{\link{LTRpred}}, \code{\link{PlotLTRAgeDistribution}}, \code{\link{PlotLTRTransposonWidthDistribution}},
#' \code{\link{PlotLTRWidthDistribution}}, \code{\link{PlotRanges}},
#' \code{\link{read.prediction}}, \code{\link{ReadLTRharvestPredictionSeqs}},
#' \code{\link{WritePredictionToFastA}}
#' @references S Steinbiss et al. Fine-grained annotation and classification of de novo predicted LTR retrotransposons. Nucl. Acids Res. (2009) 37 (21): 7002-7013.
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
    cat("\n") 
}
