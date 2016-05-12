#' @title Copy Number Quantification of predicted LTRs (solo LTR prediction)
#' @description Detect solo LTR copies and genomic locations of predicted LTR transposons using a BLAST search strategy.
#' @param LTRpred.folder file path to the \code{\link{LTRpred}} output folder.
#' @param genome file path to the reference genome in which solo LTRs shall be found (in \code{fasta} format).
#' @param ltr.similarity similarity threshold for defining LTR similarity.
#' @param scope.cutoff
#' @param perc.ident.cutoff
#' @param output file name of the BLAST output. If \code{output = NULL} (default) then the BLAST output file will be deleted after the result \code{data.frame} is returned by this function.
#' @param max.hits maximum number of hits that shall be retrieved that still fulfill the e-value criterium.
#' Default is \code{max.hits = 65000}.
#' @param eval e-value threshold for BLAST hit detection. Default is \code{eval = 1E-10}.
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

ltr.cn <- function(LTRpred.folder, 
                   genome,
                   ltr.similarity    = 70,
                   scope.cutoff      = 0.85,
                   perc.ident.cutoff = 70,
                   output            = NULL, 
                   max.hits          = 65000, 
                   eval              = 1E-10){
  
    ltr_similarity <- NULL
    
    if (is.null(output)) {
        output_3ltr <- file.path(tempdir(),"solo_ltrs_blast_output_3ltr.txt")
        output_5ltr <- file.path(tempdir(),"solo_ltrs_blast_output_5ltr.txt")
    }
    
    if (!is.null(output)) {
        output_3ltr <- file.path(tempdir(),paste0(output,"_3ltr"))
        output_5ltr <- file.path(tempdir(),paste0(output,"_5ltr"))
    }
    
    folder.name <- basename(LTRpred.folder)
    folder.name <- stringr::str_replace(folder.name,"_ltrpred","")
    
    ltrdigest.folder <- paste0(folder.name,"_ltrdigest")
    
    LTRpred.tbl <- read.ltrpred(file.path(LTRpred.folder,paste0(folder.name,"_LTRpred_DataSheet.csv")))
    
    LTR.fasta_3ltr <- file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_3ltr.fas"))
    LTR.fasta_5ltr <- file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_5ltr.fas"))
    
    LTR.filtered.fasta_3ltr <- file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_3ltr_",ltr.similarity,".fas"))
    LTR.filtered.fasta_5ltr <- file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_5ltr_",ltr.similarity,".fas"))
    
    pred2fasta(LTRpred.tbl     = dplyr::filter(LTRpred.tbl, ltr_similarity >= ltr.similarity),
               prediction.file = LTR.fasta_3ltr,
               output          = LTR.filtered.fasta_3ltr)
    
    pred2fasta(LTRpred.tbl     = dplyr::filter(LTRpred.tbl, ltr_similarity >= ltr.similarity),
               prediction.file = LTR.fasta_5ltr,
               output          = LTR.filtered.fasta_5ltr)
    
    # BLAST putative LTRs against genome file for 3 prime LTR
    system( paste0("blastn -query ",ws.wrap.path(LTR.filtered.fasta_3ltr)," -subject ",ws.wrap.path(genome),
                   " -out ", ws.wrap.path(output_3ltr) ," ","-evalue ", eval," -max_target_seqs ",max.hits," -dust no -outfmt '6 qseqid sseqid pident nident 
                   length mismatch gapopen gaps positive ppos qstart qend qlen qcovs qcovhsp qseq sstart send slen sseq evalue bitscore score'")
    )
    
    # BLAST putative LTRs against genome file for 5 prime LTR
    system( paste0("blastn -query ",ws.wrap.path(LTR.filtered.fasta_5ltr)," -subject ",ws.wrap.path(genome),
                   " -out ", ws.wrap.path(output_5ltr) ," ","-evalue ", eval," -max_target_seqs ",max.hits," -dust no -outfmt '6 qseqid sseqid pident nident 
                   length mismatch gapopen gaps positive ppos qstart qend qlen qcovs qcovhsp qseq sstart send slen sseq evalue bitscore score'")
    )
    
    # Define the column names of the BLAST output
    BLASTColNames <- c("query_id","subject_id", "perc_identity",
                       "num_ident_matches","alig_length","mismatches", "gap_openings","n_gaps","pos_match",
                       "ppos", "q_start", "q_end","q_len","qcov","qcovhsp","query_seq","s_start","s_end",
                       "s_len","subject_seq","evalue","bit_score","score_raw")
    
    # Read the BLAST output and define the columns
    BLASTOutput_3ltr <- readr::read_tsv(output_3ltr, col_names = FALSE)
    BLASTOutput_5ltr <- readr::read_tsv(output_5ltr, col_names = FALSE)
    
    colnames(BLASTOutput_3ltr) <- BLASTColNames
    colnames(BLASTOutput_5ltr) <- BLASTColNames
    
    # generate scope variable
    BLASTOutput_3ltr <- dplyr::mutate(BLASTOutput_3ltr, scope = abs(s_len - alig_length) / s_len)
    BLASTOutput_5ltr <- dplyr::mutate(BLASTOutput_5ltr, scope = abs(s_len - alig_length) / s_len)
    # filter for potentially (biologically meaningful) solo LTR hits
    BLASTOutput_3ltr <- dplyr::filter(BLASTOutput_3ltr, scope >= scope.cutoff, perc_identity >= perc.ident.cutoff) 
    BLASTOutput_5ltr <- dplyr::filter(BLASTOutput_5ltr, scope >= scope.cutoff, perc_identity >= perc.ident.cutoff) 
    # assign strand information
    BLASTOutput_3ltr <- dplyr::mutate(BLASTOutput_3ltr, strand = ifelse(s_start < s_end,"+","-")) 
    BLASTOutput_5ltr <- dplyr::mutate(BLASTOutput_5ltr, strand = ifelse(s_start < s_end,"+","-")) 
    
    # swap start and end positions for "-" strand
    BLASTOutput_3ltr <- dplyr::mutate(BLASTOutput_3ltr, s_start_new = ifelse(s_start < s_end,s_start,s_end), s_end_new = ifelse(s_start < s_end,s_end,s_start))
    BLASTOutput_5ltr <- dplyr::mutate(BLASTOutput_5ltr, s_start_new = ifelse(s_start < s_end,s_start,s_end), s_end_new = ifelse(s_start < s_end,s_end,s_start))
    BLASTOutput_3ltr <- dplyr::mutate(BLASTOutput_3ltr, s_start = s_start_new, s_end = s_end_new)
    BLASTOutput_5ltr <- dplyr::mutate(BLASTOutput_5ltr, s_start = s_start_new, s_end = s_end_new)
    
    BLASTOutput_3ltr <- dplyr::mutate(BLASTOutput_3ltr, subject_id = as.character(subject_id))
    BLASTOutput_5ltr <- dplyr::mutate(BLASTOutput_5ltr, subject_id = as.character(subject_id))
    BLASTOutput_3ltr <- dplyr::mutate(BLASTOutput_3ltr, subject_id = stringr::str_replace(subject_id,subject_id,paste0(subject_id,"_CHROMOSOME_dumped_")))
    BLASTOutput_5ltr <- dplyr::mutate(BLASTOutput_5ltr, subject_id = stringr::str_replace(subject_id,subject_id,paste0(subject_id,"_CHROMOSOME_dumped_")))
    
    # read sequences of predicted LTR transposon 
    LTR.fasta_full.te <- read.ltrpred(file.path(LTRpred.folder,paste0(folder.name,"_LTRpred_DataSheet.csv")))
    # remove mitochondria
    LTR.fasta_full.te <- dplyr::filter(LTR.fasta_full.te, !stringr::str_detect(chromosome,"mito"))
    # test whether or not 3ltr and 5 ltr loci overlap with predicted full ltr transposon locus
    full.te.chr <- names(table(LTR.fasta_full.te$chromosome))
    ltr_chr <- names(table(BLASTOutput_3ltr$subject_id))

    if (!identical(full.te.chr,ltr_chr))
        stop("Chromosome names in full LTR transposon sequence file and LTR element blast file do not match!","\n",
              "TE name = ",full.te.chr[1], " and LTR blast name = ",ltr_chr[1],". Please fix...", call. = FALSE)
    
    BLAST.ir.3ltr.list <- vector("list",length(full.te.chr))
    BLAST.ir.5ltr.list <- vector("list",length(full.te.chr))
    
    for (i in seq_len(length(full.te.chr))) {
        
        BLASTOutput_3ltr_chr <- dplyr::filter(BLASTOutput_3ltr, subject_id == full.te.chr[i])
        BLASTOutput_5ltr_chr <- dplyr::filter(BLASTOutput_5ltr, subject_id == full.te.chr[i])
        LTR.fasta_full.te_chr <- dplyr::filter(LTR.fasta_full.te, chromosome == full.te.chr[i])
        
        ir.3ltr <- IRanges::IRanges(start = BLASTOutput_3ltr_chr$s_start, end = BLASTOutput_3ltr_chr$s_end)
        ir.5ltr <- IRanges::IRanges(start = BLASTOutput_5ltr_chr$s_start, end = BLASTOutput_5ltr_chr$s_end)
        ir.full.te <- IRanges::IRanges(start = LTR.fasta_full.te_chr$start, end = LTR.fasta_full.te_chr$end)
        
        # overlap between 3ltr locus and full ltr transposon locus
        ir.3ltr_ov_ir.full.te <- IRanges::findOverlaps(ir.3ltr,ir.full.te, type = "any")
        # overlap between 5ltr locus and full ltr transposon locus
        ir.5ltr_ov_ir.full.te <- IRanges::findOverlaps(ir.5ltr,ir.full.te, type = "any")
        # overlap between 3ltr locus and 5ltr locus (do they report the same solo ltr locus?)
        ir.5ltr_ov_ir.3ltr <- IRanges::findOverlaps(ir.5ltr,ir.3ltr, type = "equal")
        # any overlap with itself: 3ltr locus with 3ltr locus
        ir.3ltr_ov_ir.3ltr <- IRanges::findOverlaps(ir.3ltr,ir.3ltr, type = "within")
        # any overlap with itself: 5ltr locus with 5ltr locus
        ir.5ltr_ov_ir.5ltr <- IRanges::findOverlaps(ir.5ltr,ir.5ltr, type = "within")
        
        nested_3ltr <- vector("list")
        for (j in names(table(ir.3ltr_ov_ir.3ltr@to))) {
            nested_3ltr_all <- which((ir.3ltr_ov_ir.3ltr@from == j) && (ir.3ltr_ov_ir.3ltr@to == j))
            nested_3ltr_subset <- ir.3ltr_ov_ir.3ltr[nested_3ltr_all, ]
            nested_3ltr[j] <- list(nested_3ltr_subset[which.max(nested_3ltr_subset[ , "bit_score"]), c("query_id","subject_id","s_start","s_end")])
        }
        
        nested_3ltr <- do.call(rbind,nested_3ltr)
        
        print(head(nested_3ltr))
        print(nrow(nested_3ltr))
        
        print(length(unique(ir.3ltr_ov_ir.3ltr@to)))
        cat("\n")
        print(length(unique(ir.5ltr_ov_ir.5ltr@to)))
        cat("\n")
        
        if ((length(ir.3ltr_ov_ir.full.te) > 0) & (length(ir.5ltr_ov_ir.3ltr) > 0) & (length(ir.5ltr_ov_ir.full.te) > 0)) {
            # exclude overlaps between 3ltr locus and full ltr transposon locus and keep equal overlaps between 3ltr locus and 5ltr locus
            BLASTOutput_3ltr_chr <- BLASTOutput_3ltr_chr[-c(ir.3ltr_ov_ir.full.te@from), ]
            # exclude overlaps between 5ltr locus and full ltr transposon locus and remove equal overlaps between 3ltr locus and 5ltr locus
            BLASTOutput_5ltr_chr <- BLASTOutput_5ltr_chr[-c(ir.5ltr_ov_ir.full.te@from,ir.5ltr_ov_ir.3ltr@from), ]
        } 
        
        if ((length(ir.3ltr_ov_ir.full.te) > 0) & (length(ir.5ltr_ov_ir.3ltr) == 0) & (length(ir.5ltr_ov_ir.full.te) > 0)) {
            # exclude overlaps between 3ltr locus and full ltr transposon locus 
            BLASTOutput_3ltr_chr <- BLASTOutput_3ltr_chr[-c(ir.3ltr_ov_ir.full.te@from), ]
            # exclude overlaps between 5ltr locus and full ltr transposon locus 
            BLASTOutput_5ltr_chr <- BLASTOutput_5ltr_chr[-c(ir.5ltr_ov_ir.full.te@from), ]
        } 
        
        if ((length(ir.3ltr_ov_ir.full.te) == 0) & (length(ir.5ltr_ov_ir.3ltr) == 0) & (length(ir.5ltr_ov_ir.full.te) > 0)) {
            # exclude overlaps between 5ltr locus and full ltr transposon locus 
            BLASTOutput_5ltr_chr <- BLASTOutput_5ltr_chr[-c(ir.5ltr_ov_ir.full.te@from), ]
        }
        
        if ((length(ir.3ltr_ov_ir.full.te) > 0) & (length(ir.5ltr_ov_ir.3ltr) == 0) & (length(ir.5ltr_ov_ir.full.te) == 0)) {
            # exclude overlaps between 3ltr locus and full ltr transposon locus and overlaps between 3ltr locus and 5ltr locus
            BLASTOutput_3ltr_chr <- BLASTOutput_3ltr_chr[-c(ir.3ltr_ov_ir.full.te@from), ]
        }
        #[ ,c("query_id","q_len","subject_id","s_start","s_end","scope","strand", "alig_length", "perc_identity", "bit_score", "evalue")]
        BLAST.ir.3ltr.list[i] <- list(dplyr::select(BLASTOutput_3ltr_chr, query_id,q_len,subject_id,s_start,s_end,scope,strand, alig_length, perc_identity, bit_score, evalue))
        BLAST.ir.5ltr.list[i] <- list(dplyr::select(BLASTOutput_5ltr_chr, query_id,q_len,subject_id,s_start,s_end,scope,strand, alig_length, perc_identity, bit_score, evalue))
    }
    
    ir.3ltr.result <- do.call(rbind,BLAST.ir.3ltr.list)
    ir.5ltr.result <- do.call(rbind,BLAST.ir.5ltr.list)
    
    res <- vector("list",2)
    res <- list(pred_3ltr = ir.3ltr.result,pred_5ltr = ir.5ltr.result)
    
    return(res)
}



