#' @title Filter the Repbase query output
#' @description Filter the output of the \code{\link{QueryRepBase}} function to quantify
#' the number of hits for each query LTR transposon (duplicates) and retain
#' only hits found in Repbase that span the annotation sequence in Repbase
#' to a certain percentage (\code{scope}).
#' @param query.output a \code{data.frame} returned by the \code{\link{QueryRepBase}} function.
#' @param scope.value a value between [0,1] qunatifying the percentage of minimum sequence similariy 
#' between the LTR transposon and the corresponding annotated sequence found in Repbase.
#' @param verbose a logical value indicating whether or not additional information shall be printed 
#' to the console while executing this function.
#' @author Hajk-Georg Drost
#' @examples 
#' \dontrun{
#' # PreProcess Repbase: A thaliana
#' # and save the output into the file "Athaliana_repbase.ref"
#' CleanRepBase(repbase.file = "athrep.ref",
#'              output.file  = "Athaliana_repbase.ref")
#'              
#' # perform blastn search against A thaliana repbase annotation
#' AthalianaRepBaseAnnotation <- QueryRepBase(ltr.seqs     = "TAIR10_chr_all-ltrdigest_complete.fas", 
#'                                            repbase.path = "Athaliana_repbase.ref", 
#'                                            cores        = 1)
#'  # filter the annotation query output                                           
#'  AthalianaAnnot.HighMatches <- FilterRepBaseQueryOutput(AthalianaRepBaseAnnotation, 
#'                                                         scope = 0.9)
#'  Ath.TE.Matches.Families <- sort(table(
#'                             unlist(lapply(stringr::str_split(
#'                             names(table(AthalianaAnnot.HighMatches$subject_id)),"_"),
#'                             function(x) paste0(x[2:3],collapse = ".")))),
#'                                         decreasing = TRUE)
#'  
#'  # visualize the hits found to have a scope of 90%
#'  barplot(Ath.TE.Matches.Families,
#'         las       = 3, 
#'         cex.names = 0.8,
#'         col       = bcolor(length(Ath.TE.Matches.Families)), 
#'         main = "RepBase Annotation: A. thaliana")
#' 
#' }
#' @return A \code{data.frame} storing the filtered output returned by \code{\link{QueryRepBase}}. 
#' @seealso \code{\link{QueryRepBase}}, \code{\link{CleanRepBase}} 


FilterRepBaseQueryOutput <- function(query.output, scope.value = 0.7, verbose = TRUE){
    
    if (!dplyr::between(scope.value,0,1))
        stop("The scope.value must be in the interval [0,1].")
    
    query_id <- bit_score <- q_len <- evalue <- scope <- NULL
  
    if (verbose){
        cat("Number of hits found by querying Repbase: ",nrow(query.output))
        cat("\n")
        cat("Number of duplicated hits found in Repbase: ",length(which(duplicated(query.output[ , "query_id"]))))
        cat("\n")
        cat("Cumulative distribution of scope: ", paste0("[ ",paste0(seq(0,1,.1) * 100, "% : "),format(quantile(unlist(query.output[ , "scope"]),probs = seq(0,1,.1)),digits = 3)," ]"))
        cat("\n")
        cat("Selecting hit with highest bit score...")
        HighestBitScoringHit <- dplyr::select(dplyr::filter(dplyr::group_by(query.output,query_id), (bit_score == max(bit_score))),query_id:q_len,evalue,bit_score,scope)
        cat("\n")
        cat("Number of remaining duplicated hits found in Repbase (due to equivalent bitscore): ",length(which(duplicated(HighestBitScoringHit[ , "query_id"]))))
        cat("\n")
        cat("Filter for hits having a scope of >= ",scope.value)
        cat("\n")
        HighestBitScoringHit <- dplyr::filter(HighestBitScoringHit, scope >= scope.value)
        cat("Cumulative distribution of scope after filtering: ", paste0("[ ",paste0(seq(0,1,.1) * 100, "% : "),format(quantile(unlist(HighestBitScoringHit[ , "scope"]),probs = seq(0,1,.1)),digits = 3)," ]"))
        cat("\n")
        cat("Number of remaining hits after filtering: ",nrow(HighestBitScoringHit))
    } else {
        
        HighestBitScoringHit <- dplyr::select(dplyr::filter(dplyr::group_by(query.output,query_id), (bit_score == max(bit_score))),query_id:q_len,evalue,bit_score,scope)
        HighestBitScoringHit <- dplyr::filter(HighestBitScoringHit, scope >= scope.value)
    }
    
    return(HighestBitScoringHit)
}






