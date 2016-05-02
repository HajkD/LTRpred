#' @title Match LTRharvest, LTRdigest, or LTRpred prediction with a given annotation file in GFF3 format
#' @description Match \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}} predictions with a given annotation file in GFF3 format
#' to check whether or not predicted LTR transposons overlap with annotated genomic features (in this case with genes).
#' @param pred.file a \code{data.frame} returned by \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}.
#' @param annotation.file imported genome annotation file from gff3 file.
#' @param strand.ori the strand orientation that shall be matched between the the prediction and the annotation.
#' Please note that some predictions do not include strand information. In case strand information should not be considered
#' please specify \code{strand.ori = NULL}.
#' @param overlap.type type of overlap between the LTR transposon prediction and the genomic feature annotated in the annotation file. 
#' Inherited from \pkg{IRanges}. Types can be \code{overlap.type = "any"}, \code{overlap.type = "start"}, \code{overlap.type = "end"},
#' \code{overlap.type = "within"}, or \code{overlap.type = "equal"}. 
#' @author Hajk-Georg Drost 
#' @examples 
#' \dontrun{
#' #' ## Hypothetical Example for Arabidopsis thaliana
#' # Generate LTR transposon prediction for A. thaliana
#'  LTRharvest("Genome/TAIR10_chr_all.fas")
#' 
#'  LTRdigest(input.gff3        = "TAIR10_chr_all/TAIR10_chr_all_Prediction.gff", 
#'            genome.file       = "Genome/TAIR10_chr_all.fas",
#'            trnas             = "araTha1-tRNAs.fa",
#'            hmms              = "hmm_*",
#'            cores             = 1)
#'
#' # Read the output of LTRdigest()
#' LTRdigest <- read.prediction(gff.file    = "TAIR10_chr_all_LTRdigestPrediction.gff",
#'                              tabout.file = "TAIR10_chr_all-ltrdigest_tabout.csv",
#'                              program     = "LTRdigest")
#' # read annotation file
#' # read the annotation file in gff3 format
#' annotation.file <- readr::read_tsv("Arabidopsis_thaliana.TAIR10.30.chr.gff3", col_names = FALSE, skip = 0, comment = "#")
#' names(annotation.file)[1:9] <- c("chromosome","source","feature","start","end","score","strand","frame","attribute")
#'                                              
#'                                        
#' # Match predicted LTR transposons with annotation file (feature = genes)        
#' Anno <- pred2annotation(
#'                pred.file = Ath.LTRdigest.prediction$ltr.retrotransposon,
#'                annotation.file       = annotation.file)
#' }
#' @export

pred2annotation <- function(pred.file,
                            annotation.file, 
                            strand.ori      = "+",
                            overlap.type    = "any" ){

    if (!is.element(strand.ori,c("+","-")))
        stop ("Please enter the correct strand.ori information! Either '+' strand or '-' strand.")
    
    chromosome <- strand <- feature <- NULL
  
    chrms.annotation <- names(table(annotation.file[ , "chromosome"]))
    chrms.pred <- names(table(pred.file[ , "chromosome"]))
    
    if (!identical(chrms.annotation,chrms.pred))
      stop ("Please check chromosome names in pred.file and annotation.file... 
            chromosome names do not seem to match.")
    
    # match LTRdigest prediction output with Annotation File in GFF3 format
    Annotation <- vector("list")
    for (i in seq_len(length(chrms.pred))){
        PutativeLTRsFiltered <- dplyr::filter(pred.file,chromosome == chrms.pred[i], strand == strand.ori) 
        GeneAnnotation <- dplyr::filter(annotation.file, (is.element(feature, c("gene","transposable_element","transposable_element_gene","transposon_fragment"))) & (chromosome == chrms.annotation[i]), strand == strand.ori)

        GeneAnnotation.bins <- IRanges::IRanges(GeneAnnotation$start, GeneAnnotation$end)
        PutativeLTRs.bins <- IRanges::IRanges(PutativeLTRsFiltered$start, PutativeLTRsFiltered$end)
        IntersectingIntervals <- IRanges::findOverlaps(PutativeLTRs.bins,GeneAnnotation.bins, type = overlap.type)
        Annotation[i] <- list(cbind(PutativeLTRsFiltered[IntersectingIntervals@queryHits, ], GeneAnnotation[IntersectingIntervals@subjectHits, ]))
    }
    AnnotationResult <- do.call(rbind,Annotation)
    #AnnotationResult <- AnnotationResult[order(AnnotationResult[ , "ltr_similarity"], decreasing = TRUE), ]
    
    return (AnnotationResult)
}

