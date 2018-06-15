#' @title Estimate retrotransposon insertion age in Mya based on 5' and 3' LTR homology
#' @description This function implements diverse metrics to roughly estimate
#' the insertion age in Mya based on 5' and 3' LTR homology.
#' @param ltr_seqs_3_prime file path to a fasta file storing the sequences of the respective 3' LTR.
#' @param ltr_seqs_5_prime file path to a fasta file storing the sequences of the respective 5' LTR.
#' @param model a model as specified in \code{\link[ape]{dist.dna}}: a character string specifying the evolutionary model to be used; must be one of "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin".
#' @param mutation_rate a mutation rate per site per year. For retrotransposons the default is \code{mutation_rate = 1.3 * 10E-8} (Wicker and Keller, 2007).
#' @author Hajk-Georg Drost
#' @export

ltr_age_estimation <-
  function(ltr_seqs_3_prime,
           ltr_seqs_5_prime,
           model = "K80",
           mutation_rate = 1.3 * 10E-8) {
    
    
    # import 3' and 5' LTR sequences
    ltr_3_prime_seqs <- Biostrings::readDNAStringSet(ltr_seqs_3_prime)
    ltr_5_prime_seqs <- Biostrings::readDNAStringSet(ltr_seqs_5_prime)
    
    if (length(ltr_3_prime_seqs) != length(ltr_5_prime_seqs))
      stop("Something went wrong! The files '",ltr_seqs_3_prime,"' and '",ltr_seqs_5_prime,"' don't contain the same number of LTR sequences.", call. = FALSE)
    
    # compute global pairwise alignments between 3' and 5' LTR sequences
    Alignments <-
      Biostrings::pairwiseAlignment(ltr_3_prime_seqs, ltr_5_prime_seqs)
    
    names_3_prime_seqs <- names(ltr_3_prime_seqs)
    names_5_prime_seqs <- names(ltr_5_prime_seqs)
    
    # retrieve 3' LTR alignment sequences
    ltr_3_prime_seqs_with_gaps <-
      Biostrings::DNAStringSet(Alignments@pattern)
    names(ltr_3_prime_seqs_with_gaps) <- names(ltr_3_prime_seqs)
    # save 3' LTR alignment sequences for import into the ape package
    Biostrings::writeXStringSet(
      ltr_3_prime_seqs_with_gaps,
      file.path(tempdir(), "ltr_3_prime_seqs_with_gaps.fasta"), 
    )
    # import 3' LTR alignment sequences to ape
    ltr_3_prime_seqs_with_gaps_ape <-
      ape::read.dna(file.path(tempdir(), "ltr_3_prime_seqs_with_gaps.fasta"),
                    format = "fasta")
    
    # retrieve 5' LTR alignment sequences
    ltr_5_prime_seqs_with_gaps <-
      Biostrings::DNAStringSet(Alignments@subject)
    names(ltr_5_prime_seqs_with_gaps) <- names(ltr_5_prime_seqs)
    # save 5' LTR alignment sequences for import into the ape package
    Biostrings::writeXStringSet(
      ltr_5_prime_seqs_with_gaps,
      file.path(tempdir(), "ltr_5_prime_seqs_with_gaps.fasta")
    )
    
    # import 5' LTR alignment sequences to ape
    ltr_5_prime_seqs_with_gaps_ape <-
      ape::read.dna(file.path(tempdir(), "ltr_5_prime_seqs_with_gaps.fasta"),
                    format = "fasta")
    
    print(ltr_5_prime_seqs_with_gaps_ape)
    
    dist_vals <- vector("numeric", length(ltr_5_prime_seqs_with_gaps_ape))
    indel_blocks <- vector("numeric", length(ltr_5_prime_seqs_with_gaps_ape))
    indel <- vector("numeric", length(ltr_5_prime_seqs_with_gaps_ape))
    
    for (i in seq_len(length(ltr_5_prime_seqs_with_gaps_ape))) {
      # print(ltr_5_prime_seqs_with_gaps_ape[[i]])
      seq_5_prime <- Biostrings::DNAStringSet(ltr_5_prime_seqs_with_gaps_ape[[i]])
      names(seq_5_prime) <- names_5_prime_seqs[i]
      
      Biostrings::writeXStringSet(
        Biostrings::DNAStringSet(seq_5_prime),
        file.path(tempdir(), paste0("seq_",i, ".fasta"))
      )
      
      seq_3_prime <- Biostrings::DNAStringSet(ltr_3_prime_seqs_with_gaps_ape[[i]])
      names(seq_3_prime) <- names_3_prime_seqs[i]
      
      Biostrings::writeXStringSet(
        seq_3_prime,
        file.path(tempdir(), paste0("seq_",i, ".fasta")),
        append = TRUE
      )
      
      seq_i <-
        ape::read.dna(paste0("seq_",i, ".fasta"), format = "fasta")
      dist_vals[i] <- ape::dist.dna(seq_i, model = model, variance = TRUE)
      indel_blocks[i] <- ape::dist.dna(seq_i, model = "indelblock")
      indel[i] <- ape::dist.dna(seq_i, model = "indel")
    }
    
    res <- tibble::tibble(
      ltr_name = names(ltr_3_prime_seqs),
      ltr_age_mya = (dist_vals / (2 * mutation_rate)) / 1E6,
      ltr_aln_score = Alignments@score,
      ltr_indel_blocks = indel_blocks,
      ltr_indels = indel
    )
    
    res <- dplyr::mutate(res, ltr_age_naiv_mya = res$ltr_indels * res$mutation_rate * 1E6)
    
    return(res)
  }
