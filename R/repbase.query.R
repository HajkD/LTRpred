#' @title Query the RepBase to annotate putative LTRs
#' @description Validate or annotate putative LTR
#' transposons that have been predicted using LTRharvest or LTRdigest.
#' @param seq.file file path to the putative LTR transposon sequences in \code{fasta} format.
#' @param repbase.path file path to the RepBase file in \code{fasta} format.
#' @param output file name of the BLAST output.
#' @param max.hits maximum number of hits that shall be retrieved that still fulfill the e-value criterium.
#' Default is \code{max.hits = 5000}.
#' @param eval e-value threshold for BLAST hit detection. Default is \code{eval = 1E-30}.
#' @param cores number of cores to use to perform parallel computations.
#' @author Hajk-Georg Drost
#' @details 
#' The RepBase database provides a collection of curated transposable element annotations.
#' 
#' This function allows users to validate or annotate putative LTR
#' transposons that have been predicted using LTRharvest or LTRdigest by blasting predicted LTR transposons 
#' to transposons known (annotated) in other species (e.g. such as Arabidopsis thaliana).
#' 
#' Internally, this function performs a \code{blastn} search of the putative LTR transposons predicted
#' by LTRharvest or LTRdigest against the Repbase fasta file that is specified by the user.
#' 
#' For this purpose it is required that the user has a working version of BLAST+ running on his or her machine.
#' @examples 
#' \dontrun{
#' # Example annotation run against the A thaliana RepBase using 4 cores
#' q <- repbase.query(seq.file     = "path/to/LTRtransposonSeqs.fasta",
#'                   repbase.path = "path/to/Athaliana_repbase.ref",
#'                   cores        = 4)
#'                  
#' Annot <- dplyr::select(dplyr::filter(dplyr::group_by(q,query_id), 
#'                                     (bit_score == max(bit_score))),
#'                                      query_id:q_len,evalue,bit_score,scope)
#' # select only hits with a scope > 0.1
#' Annot.HighMatches <- dplyr::filter(Annot, scope >= 0.1)
#' # Annotate the proportion of hits
#' barplot(sort(table(unlist(lapply(stringr::str_split(
#'         names(table(Annot.HighMatches$subject_id)),"_"), 
#'         function(x) x[2]))), decreasing = TRUE))
#' }
#' @references 
#' http://www.girinst.org/repbase/
#' 
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
      
repbase.query <- function(seq.file, 
                         repbase.path, 
                         output   = "RepbaseOutput.txt", 
                         max.hits = 5000, 
                         eval     = 1E-30, 
                         cores    = 1){
    
    s_len <- alig_length <- NULL
    # Create a blast-able database of repbase and perform blastn search of 
    # LTR transposons against this blast formatted repbase file
    system(paste0("makeblastdb -in ",repbase.path," -parse_seqids -input_type fasta -dbtype nucl"))
    system( paste0("blastn -db ",repbase.path," -query ",seq.file,
                   " -out ", output ," ","-evalue ", eval," -max_target_seqs ",max.hits," -num_threads ",cores,
                   " -dust no -outfmt '6 qseqid sseqid pident nident 
                   length mismatch gapopen gaps positive ppos qstart qend qlen qcovs qcovhsp qseq sstart send slen sseq evalue bitscore score'")
    )
    
    # Define the column names of the BLAST output
    BLASTColNames <- c("query_id","subject_id", "perc_identity",
                       "num_ident_matches","alig_length","mismatches", "gap_openings","n_gaps","pos_match",
                       "ppos", "q_start", "q_end","q_len","qcov","qcovhsp","query_seq","s_start","s_end",
                       "s_len","subject_seq","evalue","bit_score","score_raw")
    # Read the BLAST output and define the columns
    BLASTRepBaseOutput <- readr::read_tsv(output, col_names = FALSE)
    colnames(BLASTRepBaseOutput) <- BLASTColNames
    # Add an additional variable named 'scope' to the BLAST output
    # Scope is defined by the formula: 1 - (abs(s_len - alig_length) / s_len)
    # and is computed for each hit of the BLAST output
    # This scope measure allows to quantify the percentage of sequence similarity
    # over the entire query and subject sequence length.
    # A scope of 0.5 can be interpreted as a 50% coverage (in terms of sequence similarity)
    # found between the query and subject DNA string
    BLASTRepBaseOutput <- dplyr::mutate(BLASTRepBaseOutput, scope = 1 - (abs(s_len - alig_length) / s_len))
    
    return(BLASTRepBaseOutput)
}

