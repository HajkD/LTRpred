#' @title Low level function to detect motifs in strings
#' @description Find a specific motif or a sequence of motifs within
#' genomic sequences.
#' @param seq.file path to the genomic sequecne file of interest (e.g. LTR TE seqs predicted by
#' \code{\link{LTRpred}}).
#' @param motif a character string or vector of strings which shall be counted within each sequence.
#' @param as.ratio shall count values be returned as asbolute frequency (count value) or as relative frequency (percentage).
#' @author Hajk-Georg Drost
#' @examples 
#' # find number of "CG" motifs in predicted LTR transposons
#' motif.count(seq.file = system.file("LTRseqs.fas",package = "LTRpred"), 
#'             motif    = "CG")
#'             
#' # find number of "CG" motifs in predicted LTR transposons: rel. frequency
#' motif.count(seq.file = system.file("LTRseqs.fas",package = "LTRpred"), 
#'             motif    = "CG",
#'             as.ratio = TRUE)             
#' @export
       
motif.count <- function(seq.file, motif, as.ratio = FALSE){
   
    # read sequence
    seqs <- Biostrings::readDNAStringSet(seq.file)
    
    if (as.ratio) {
        if (motif == "CH") {
            seq_lengths <- vector("numeric", length(seqs))
            res.CA <- vector("numeric", length(seqs))
            res.CC <- vector("numeric", length(seqs))
            res.CT <- vector("numeric", length(seqs))
            res <- vector("numeric", length(seqs))
            
            seq_lengths <-  Biostrings::nchar(seqs)
            res.CA <-
                as.numeric(Biostrings::vcountPattern("CA", seqs)) / seq_lengths
            res.CC <-
                as.numeric(Biostrings::vcountPattern("CC", seqs)) / seq_lengths
            res.CT <-
                as.numeric(Biostrings::vcountPattern("CT", seqs)) / seq_lengths
            res <- res.CA + res.CC + res.CT
        }
        
        if (motif == "CHG") {
            seq_lengths <- vector("numeric", length(seqs))
            seq_lengths <-  Biostrings::nchar(seqs)
            res.CAG <- vector("numeric", length(seqs))
            res.CCG <- vector("numeric", length(seqs))
            res.CTG <- vector("numeric", length(seqs))
            res <- vector("numeric", length(seqs))
            
            res.CAG <-
                as.numeric(Biostrings::vcountPattern("CAG", seqs)) / seq_lengths
            res.CCG <-
                as.numeric(Biostrings::vcountPattern("CCG", seqs)) / seq_lengths
            res.CTG <-
                as.numeric(Biostrings::vcountPattern("CTG", seqs)) / seq_lengths
            res <- res.CAG + res.CCG + res.CTG
        }
        
        if (motif == "CHH") {
            seq_lengths <- vector("numeric", length(seqs))
            seq_lengths <-  Biostrings::nchar(seqs)
            res <- vector("numeric", length(seqs))
            
            # retrieve all combinations of CHH
            CHH.comb <-
                apply(expand.grid(list(
                    "C", c("A", "C", "T"), c("A", "C", "T")
                )), 1, stringr::str_c, collapse = "")
            rel.count <-
                function(x)
                    as.numeric(Biostrings::vcountPattern(x, seqs))
            res <-
                rowSums(sapply(CHH.comb, rel.count)) / seq_lengths
        }
        
        if (!is.element(motif, c("CH", "CHG", "CHH"))) {
            seq_lengths <- vector("numeric", length(seqs))
            seq_lengths <-  Biostrings::nchar(seqs)
            res <- vector("numeric", length(seqs))
            
            res <-
                as.numeric(Biostrings::vcountPattern(motif, seqs)) / seq_lengths
        }
        
        names(res) <- seqs@ranges@NAMES
    }
    
    if (!as.ratio) {
        if (motif == "CH") {
            seq_lengths <- vector("numeric", length(seqs))
            res.CA <- vector("numeric", length(seqs))
            res.CC <- vector("numeric", length(seqs))
            res.CT <- vector("numeric", length(seqs))
            res <- vector("numeric", length(seqs))
            
            seq_lengths <-  as.numeric(Biostrings::nchar(seqs))
            res.CA <-
                as.numeric(Biostrings::vcountPattern("CA", seqs))
            res.CC <-
                as.numeric(Biostrings::vcountPattern("CC", seqs))
            res.CT <-
                as.numeric(Biostrings::vcountPattern("CT", seqs))
            res <- res.CA + res.CC + res.CT
        }
        
        if (motif == "CHG") {
            res.CAG <- vector("numeric", length(seqs))
            res.CCG <- vector("numeric", length(seqs))
            res.CTG <- vector("numeric", length(seqs))
            res <- vector("numeric", length(seqs))
            
            res.CAG <- as.numeric(Biostrings::vcountPattern("CAG", seqs))
            res.CCG <- as.numeric(Biostrings::vcountPattern("CCG", seqs))
            res.CTG <- as.numeric(Biostrings::vcountPattern("CTG", seqs))
            res <- res.CAG + res.CCG + res.CTG
        }
        
        if (motif == "CHH") {
            res <- vector("numeric", length(seqs))
            
            # retrieve all combinations of CHH
            CHH.comb <-
                apply(expand.grid(list(
                    "C", c("A", "C", "T"), c("A", "C", "T")
                )), 1, stringr::str_c, collapse = "")
            rel.count <-
                function(x)
                    as.numeric(Biostrings::vcountPattern(x, seqs))
            res <- rowSums(sapply(CHH.comb, rel.count))
        }
        
        if (!is.element(motif, c("CH", "CHG", "CHH"))) {
            res <- vector("numeric", length(seqs))
            res <- as.numeric(Biostrings::vcountPattern(motif, seqs))
        }
    }
  
  return(res)
}






