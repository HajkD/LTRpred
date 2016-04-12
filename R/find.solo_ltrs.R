#' @title Detect solo LTR copies of predicted LTR transposons
#' @description Detect solo LTR copies and genomic locations of predicted LTR transposons using a BLAST search strategy.
#' @param LTRpred.folder file path to the \code{\link{LTRpred}} output folder.
#' @param genome file path to the reference genome in which solo LTRs shall be found (in \code{fasta} format).
#' @param ltr.similarity similarity threshold for defining LTR similarity.
#' @param output file name of the BLAST output. If \code{output = NULL} (default) then the BLAST output file will be deleted after the result \code{data.frame} is returned by this function.
#' @param max.hits maximum number of hits that shall be retrieved that still fulfill the e-value criterium.
#' Default is \code{max.hits = 5000}.
#' @param eval e-value threshold for BLAST hit detection. Default is \code{eval = 1E-5}.
#' @author Hajk-Georg Drost
#' @details 
#' The
#' @examples 
#' \dontrun{
#' 
#' }
#' @references  
#' Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
#' 
#' Gish, W. & States, D.J. (1993) "Identification of protein coding regions by database similarity search." Nature Genet. 3:266-272.
#'
#' Madden, T.L., Tatusov, R.L. & Zhang, J. (1996) "Applications of network BLAST server" Meth. Enzymol. 266:131-141.
#'
#' Altschul, S.F., Madden, T.L., Schaeffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402.
#'
#' Zhang Z., Schwartz S., Wagner L., & Miller W. (2000), "A greedy algorithm for aligning DNA sequences" J Comput Biol 2000; 7(1-2):203-14.
#' @export

find.solo_ltrs <- function(LTRpred.folder, 
                          genome,
                          ltr.similarity = 98,
                          output   = NULL, 
                          max.hits = 5000, 
                          eval     = 1E-5){
  
  ltr_similarity <- NULL
  
  if (is.null(output))
    output <- file.path(tempdir(),"solo_ltrs_blast_output.txt")
  
  folder.name <- basename(LTRpred.folder)
  folder.name <- stringr::str_replace(folder.name,"_ltrpred","")
  
  ltrdigest.folder <- paste0(folder.name,"_ltrdigest")
  
  LTRpred.tbl <- read.ltrpred(file.path(LTRpred.folder,paste0(folder.name,"_LTRpred_DataSheet.csv")))
  
  LTR.fasta <- file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_3ltr.fas"))
  
  LTR.filtered.fasta <- file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_3ltr_",ltr.similarity,".fas"))
  
  pred2fasta(LTRpred.tbl     = dplyr::filter(LTRpred.tbl, ltr_similarity >= ltr.similarity),
             prediction.file = LTR.fasta,
             output          = LTR.filtered.fasta)
  
  # BLAST putative LTRs against genome file
  system( paste0("blastn -query ",ws.wrap.path(LTR.filtered.fasta)," -subject ",ws.wrap.path(genome),
                 " -out ", output ," ","-evalue ", eval," -max_target_seqs ",max.hits," -dust no -outfmt '6 qseqid sseqid pident nident 
                 length mismatch gapopen gaps positive ppos qstart qend qlen qcovs qcovhsp qseq sstart send slen sseq evalue bitscore score'")
  )
  
  # Define the column names of the BLAST output
  BLASTColNames <- c("query_id","subject_id", "perc_identity",
                     "num_ident_matches","alig_length","mismatches", "gap_openings","n_gaps","pos_match",
                     "ppos", "q_start", "q_end","q_len","qcov","qcovhsp","query_seq","s_start","s_end",
                     "s_len","subject_seq","evalue","bit_score","score_raw")
  # Read the BLAST output and define the columns
  BLASTOutput <- readr::read_tsv(output, col_names = FALSE)
  colnames(BLASTOutput) <- BLASTColNames
  
  
  return(BLASTOutput)
}