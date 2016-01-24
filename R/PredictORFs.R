#' @title Open Reading Frame Prediction
#' @description This function provides a wrapper to the \code{USEARCH}
#' \code{fastx_findorfs} command line tool to predict ORFs in
#' a given input fasta file and read the output as \code{\link{data.frame}} object.
#' @param input.file
#' @param orf.style
#' @param min.codons
#' @param trans.seqs
#' @param output
#' @author Hajk-Georg Drost
#' @examples 
#' \dontrun{
#' 
#' }
#' 
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
                    " -ntout ",file.path(output,paste0(basename(input.file),"_nt.fsa")),
                    " -aaout ",file.path(output,paste0(basename(input.file),"_aa.fsa")),
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



