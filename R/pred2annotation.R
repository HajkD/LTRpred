#' @title Match LTRharvest, LTRdigest, or LTRpred prediction with a given annotation file in GFF3 format
#' @description Match \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}} predictions with a given annotation file in GFF3 format
#' to check whether or not predicted LTR transposons overlap with annotated genomic features (in this case with genes).
#' @param pred.file a \code{data.frame} returned by \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, or \code{\link{LTRpred}}.
#' @param annotation.file an genome annotation file in gff3 format.
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
#'                                             
#' # Match predicted LTR transposons with annotation file (feature = genes)        
#' Anno <- pred2annotation(
#'                pred.file = Ath.LTRdigest.prediction$ltr.retrotransposon,
#'                annotation.file       = "Arabidopsis_thaliana.TAIR10.30.chr.gff3")
#' }
#' @export

pred2annotation <- function(pred.file,
                            annotation.file, 
                            strand.ori      = "+",
                            overlap.type    = "any" ){

    if (!is.element(strand.ori,c("+","-")))
        stop ("Please enter the correct strand.ori information! Either '+' strand or '-' strand.")
    
    chromosome <- strand <- feature <- seqname <- NULL
  
    # read the annotation file in gff3 format
    AnnotationFileGFF3 <- readr::read_tsv(annotation.file, col_names = FALSE, skip = 0, comment = "#")
    names(AnnotationFileGFF3)[1:9] <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
    
    # match LTRdigest prediction output with Annotation File in GFF3 format
    Annotation <- vector("list")
    for (i in as.numeric(names(table(pred.file[ , "chromosome"])))){
        PutativeLTRsFiltered <- dplyr::filter(pred.file,chromosome == i, strand == strand.ori) 
        GeneAnnotation <- dplyr::filter(AnnotationFileGFF3, (feature == "gene") & (seqname == i), strand == strand.ori)

        GeneAnnotation.bins <- IRanges::IRanges(GeneAnnotation$start, GeneAnnotation$end)
        PutativeLTRs.bins <- IRanges::IRanges(PutativeLTRsFiltered$start, PutativeLTRsFiltered$end)
        IntersectingIntervals <- IRanges::findOverlaps(PutativeLTRs.bins,GeneAnnotation.bins, type = overlap.type)
        Annotation[i] <- list(cbind(PutativeLTRsFiltered[IntersectingIntervals@queryHits, ], GeneAnnotation[IntersectingIntervals@subjectHits, ]))
    }
    AnnotationResult <- do.call(rbind,Annotation)
    AnnotationResult <- AnnotationResult[order(AnnotationResult[ , "ltr_similarity"], decreasing = TRUE), ]
    
    return (AnnotationResult)
}

