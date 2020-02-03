#' @title Estimate retrotransposon insertion age in Mya based on 5 prime and 3 prime LTR sequence homology
#' @description This function implements diverse metrics to roughly estimate
#' the insertion age in Mya based on 5 prime and 3 prime LTR sequence homology.
#' @param pred a prediction file generated with \code{\link{LTRpred}}.
#' @param ltr_seqs_3_prime file path to a fasta file storing the sequences of the respective 3 prime LTR (e.g. as annotatted by \code{\link{LTRpred}}).
#' @param ltr_seqs_5_prime file path to a fasta file storing the sequences of the respective 5 prime LTR (e.g. as annotatted by \code{\link{LTRpred}}).
#' @param model a model as specified in \code{\link[ape]{dist.dna}}: a character string specifying the evolutionary model to be used - must be one of
#'  \itemize{
#' \item  \code{K80} (the default)
#' \item \code{raw}
#' \item  \code{N}
#' \item  \code{TS}
#' \item  \code{TV}
#' \item  \code{JC69}
#' \item  \code{F81} 
#' \item \code{K81}
#' \item \code{F84}
#' \item \code{BH87}
#' \item \code{T92}
#' \item \code{TN93}
#' \item \code{GG95}
#' \item \code{logdet}
#' \item \code{paralin}
#' }
#' @param mutation_rate a mutation rate per site per year. For retrotransposons the default is \eqn{mutation_rate = 1.3 * 10E-8} (Wicker and Keller, 2007).
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' ltr_pred <- LTRpred::read.ltrpred(
#'           system.file("Hsapiens_ChrY_LTRpred_DataSheet.tsv", 
#'           package = "LTRpred"))
#' # define file path to fasta file storing 3 prime LTR sequences
#' ltr_seqs_3_prime <- system.file("Hsapiens_ChrY-ltrdigest_3ltr.fas", package = "LTRpred")
#' ltr_seqs_5_prime <- system.file("Hsapiens_ChrY-ltrdigest_5ltr.fas", package = "LTRpred")
#' # estimate insertion age based on 3 prime and 5 prime LTR homology using the K80 model
#' Hsapiens_ltr_age <- LTRpred::ltr_age_estimation(ltr_pred, ltr_seqs_3_prime, ltr_seqs_5_prime)
#' # look at results
#' Hsapiens_ltr_age
#' }
#' @export

ltr_age_estimation <-
  function(pred,
           ltr_seqs_3_prime,
           ltr_seqs_5_prime,
           model = "K80",
           mutation_rate = 1.3 * 10E-8) {
    
    message("Starting retrotransposon evolutionary age estimation by comparing the 3' and 5' LTRs using the molecular evolution model '", model, "' and the mutation rate '", mutation_rate, "' (please make sure the mutation rate can be assumed for your species of interest!) for ", nrow(pred), " predicted elements ...")
    message("\n")
    message("Please be aware that evolutionary age estimation based on 3' and 5' LTR comparisons are only very rough time estimates and don't take reverse-transcription mediated retrotransposon recombination between family members of retroelements into account! Please consult Sanchez et al., 2017 Nature Communications and Drost & Sanchez, 2019 Genome Biology and Evolution for more details on retrotransposon recombination.")
    
    if (!file.exists(ltr_seqs_3_prime))
      stop("The file '", ltr_seqs_3_prime, "' does not seem to exist. Please provide a valid file path for argument 'ltr_seqs_3_prime'.", call. = FALSE)
    
    if (!file.exists(ltr_seqs_5_prime))
      stop("The file '", ltr_seqs_5_prime, "' does not seem to exist. Please provide a valid file path for argument 'ltr_seqs_5_prime'.", call. = FALSE)
    
    new_ltr_seqs_3_prime <- file.path(tempdir(), basename(ltr_seqs_3_prime))
    new_ltr_seqs_5_prime <- file.path(tempdir(), basename(ltr_seqs_5_prime))
    pred2fasta(pred, ltr_seqs_3_prime, new_ltr_seqs_3_prime)
    pred2fasta(pred, ltr_seqs_5_prime, new_ltr_seqs_5_prime)
    
    # import 3' and 5' LTR sequences
    ltr_3_prime_seqs <- Biostrings::readDNAStringSet(new_ltr_seqs_3_prime)
    ltr_5_prime_seqs <- Biostrings::readDNAStringSet(new_ltr_seqs_5_prime)
    
    if (length(ltr_3_prime_seqs) != length(ltr_5_prime_seqs))
      stop("Something went wrong! The files '",ltr_seqs_3_prime,"' and '",ltr_seqs_5_prime,"' don't contain the same number of LTR sequences.", call. = FALSE)
    
    # compute global pairwise alignments between 3' and 5' LTR sequences
    Alignments <-
      Biostrings::pairwiseAlignment(ltr_3_prime_seqs, ltr_5_prime_seqs)
    
    # retrieve 3' LTR alignment sequences
    ltr_3_prime_seqs_with_gaps <-
      Biostrings::BStringSet(Alignments@pattern)
    names(ltr_3_prime_seqs_with_gaps) <- names(ltr_3_prime_seqs)
    
    # save 3' LTR alignment sequences for import into the ape package
    Biostrings::writeXStringSet(
      ltr_3_prime_seqs_with_gaps,
      file.path(tempdir(), paste0(basename(ltr_seqs_3_prime), "_ltr_3_prime_seqs_with_gaps.fasta")) 
    )
    # import 3' LTR alignment sequences to ape
    ltr_3_prime_seqs_with_gaps_ape <-
      ape::read.dna(file.path(tempdir(), paste0(basename(ltr_seqs_3_prime), "_ltr_3_prime_seqs_with_gaps.fasta")),
                    format = "fasta")
    
    # retrieve 5' LTR alignment sequences
    ltr_5_prime_seqs_with_gaps <-
      Biostrings::BStringSet(Alignments@subject)
    names(ltr_5_prime_seqs_with_gaps) <- names(ltr_5_prime_seqs)
    
    # save 5' LTR alignment sequences for import into the ape package
    Biostrings::writeXStringSet(
      ltr_5_prime_seqs_with_gaps,
      file.path(tempdir(), paste0(basename(ltr_seqs_5_prime), "_ltr_5_prime_seqs_with_gaps.fasta"))
    )
    
    # import 5' LTR alignment sequences to ape
    ltr_5_prime_seqs_with_gaps_ape <-
      ape::read.dna(file.path(tempdir(), paste0(basename(ltr_seqs_5_prime), "_ltr_5_prime_seqs_with_gaps.fasta")),
                    format = "fasta")
    
    dist_vals <- vector("numeric", length(ltr_5_prime_seqs_with_gaps_ape))
    indel_blocks <- vector("numeric", length(ltr_5_prime_seqs_with_gaps_ape))
    indel <- vector("numeric", length(ltr_5_prime_seqs_with_gaps_ape))
    
    for (i in seq_len(length(ltr_5_prime_seqs_with_gaps_ape))) {
      # seq_5_prime <- Biostrings::DNAStringSet(ltr_5_prime_seqs_with_gaps_ape[[i]])
      # names(seq_5_prime) <- names_5_prime_seqs[i]
      
      ape::write.FASTA(
        ltr_5_prime_seqs_with_gaps_ape[i],
        file.path(tempdir(), paste0(basename(ltr_seqs_5_prime), "_seq_",i, ".fasta"))
      )
      
      # seq_3_prime <- Biostrings::DNAStringSet(ltr_3_prime_seqs_with_gaps_ape[[i]])
      # names(seq_3_prime) <- names_3_prime_seqs[i]
      
      ape::write.FASTA(
        ltr_3_prime_seqs_with_gaps_ape[i],
        file.path(tempdir(), paste0(basename(ltr_seqs_5_prime), "_seq_",i, ".fasta")),
        append = TRUE
      )
      
      seq_i <-
        ape::read.dna(file.path(tempdir(), paste0(basename(ltr_seqs_5_prime), "_seq_",i, ".fasta")), format = "fasta")
      dist_vals[i] <- ape::dist.dna(seq_i, model = model, variance = TRUE)
      indel_blocks[i] <- ape::dist.dna(seq_i, model = "indelblock")
      indel[i] <- ape::dist.dna(seq_i, model = "indel")
    }
    
    res <- tibble::tibble(
      orf.id = pred$orf.id,
      ltr_name = names(ltr_3_prime_seqs),
      ltr_age_mya = (dist_vals / (2 * mutation_rate)) / 1E6,
      ltr_evo_distance = dist_vals,
      ltr_aln_score = Alignments@score,
      ltr_indel_blocks = indel_blocks,
      ltr_indels = indel
    )
    
    res <-
      dplyr::mutate(res, ltr_age_naive_mya = res$ltr_indels * mutation_rate * 1E6)
    
    return(res)
  }
