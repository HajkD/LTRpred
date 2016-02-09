#' @title Open Reading Frame Prediction
#' @description This function provides a wrapper to the \code{USEARCH}
#' \code{fastx_findorfs} command line tool to predict ORFs in
#' a given input fasta file and read the output as \code{\link{data.frame}} object.
#' @param input.file a fasta file storing the sequences for which open reading frames shall be predicted.
#' @param orf.style type of predicting open reading frames (see documentation of USEARCH).
#' @param min.codons minimum number of codons in the predicted open reading frame.
#' @param trans.seqs logical value indicating wheter or not predicted open reading frames
#' shall be translated and the corresponding protein sequences stored in the output folder. 
#' @param output path to the folder in which predicted open reading frame sequences shall be stored.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{LTRharvest}}, \code{\link{LTRdigest}}, \code{\link{LTRpred}}
#' @references 
#' Robert Edgar. Search and clustering orders of magnitude faster than BLAST. Bioinformatics (2010) 26 (19): 2460-2461.
#' @export

PredictORFs <- function(input.file, 
                        orf.style  = 7, 
                        min.codons = 200, 
                        trans.seqs = FALSE,
                        output     = NULL){
  
  
  if (!trans.seqs){
    
    if (is.null(output)){
      system(paste0("usearch -fastx_findorfs ",input.file,
                    " -ntout ",paste0(basename(input.file),"_nt.fsa"),
                    " -orfstyle ",orf.style," -mincodons ",min.codons))
    } else {
      
      system(paste0("usearch -fastx_findorfs ",input.file,
                    " -ntout ",file.path(output,paste0(basename(input.file),"_nt.fsa")),
                    " -orfstyle ",orf.style," -mincodons ",min.codons))
    }
    
  } else {
    
    if (is.null(output)){
      system(paste0("usearch -fastx_findorfs ",input.file,
                    " -ntout ",paste0(basename(input.file),"_nt.fsa"),
                    " -aaout ",paste0(basename(input.file),"_aa.fsa"),
                    " -orfstyle ",orf.style," -mincodons ",min.codons))
    } else {
      system(paste0("usearch -fastx_findorfs ",input.file,
                    " -ntout ",file.path(output,paste0(basename(input.file),"_orfs_nt.fsa")),
                    " -aaout ",file.path(output,paste0(basename(input.file),"_orfs_aa.fsa")),
                    " -orfstyle ",orf.style," -mincodons ",min.codons))
    }
  }
  
  if (!is.null(output)){
    ORFCount.df <- read.findorfs.output(file.path(output,paste0(basename(input.file),"_nt.fsa")))
  } else {
    ORFCount.df <- read.findorfs.output(paste0(basename(input.file),"_nt.fsa"))
  }
  
  return (ORFCount.df)
}

