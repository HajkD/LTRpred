#' @title Generating genome summary files for \code{LTRpred.meta} results
#' @details ...
#' @param genome.folder
#' @param ltrpred.meta.folder
#' @param file.name
#' @param sim
#' @param cut.range
#' @param quality.filter
#' @param n.orfs
#' @param strategy
#' @author Hajk-Georg Drost
#' @export

genome.summary <- function(genome.folder,
                           ltrpred.meta.folder,
                           file.name,
                           sim                 = 70,
                           cut.range           = 2,
                           quality.filter      = TRUE,
                           n.orfs              = 0,
                           strategy            = "default") {
    
    
    
    if (!dir.exists(genome.folder))
        stop("Genome folder '", genome.folder, "' does not exist. Please provide a valid path to the folder storing your genome assembly files.", call. = FALSE)
    
    if (!dir.exists(ltrpred.meta.folder))
        stop("The folder '", ltrpred.meta.folder, "' that is supposed to store the results of the LTRpred.meta() run does not exist! Please provide a valid path.", call. = FALSE)
    
    if (length(rev(seq(100, sim,-cut.range))) == 1L)
        stop(
            "Please specify a 'cut.range' value that is compatible with the 'sim' threshold. ",
            "The input 'cut.range = ",
            cut.range,
            "', whereas the similarity bin is between: [",
            sim,
            ",100].",
            call. = FALSE
        )
    
    assembly_files <- list.files(genome.folder)
    ltrpred_meta_results <- list.files(ltrpred.meta.folder)
    
    if (any(stringr::str_detect(assembly_files, "doc_")) | any(stringr::str_detect(assembly_files, "md5checksum")))
        stop("It seems that you downloaded genome assemblies using the biomartr package. Many thanks for using biomartr :-) .",
             " Unfortunately, your genome assembly folder still stores the 'doc_' or 'md5checksum' files. Please remove them from the folder so that ",
             "this function can run properly (only fasta files should be retained in the genome.file folder).", call. = FALSE)
    
    
    assembly_files_chop <- str_chop_vec(assembly_files, pattern = "[.]")
    ltrpred_meta_results_chop <- str_chop_vec(ltrpred_meta_results, pattern = "_")
    
    if (length(dplyr::setdiff(assembly_files_chop, ltrpred_meta_results_chop)) > 0) {
        message("\n")
        message("Please make sure that genome assembly file names in '", genome.folder, 
                "' match with the LTRpred.meta() results folder names in '",ltrpred.meta.folder,"'.",
                " E.g. Hsapiens.fa (assembly name in genome.folder) and Hsapiens_ltrpred (LTRpred.meta() result in ltrpred.meta.folder).")
        message("\n")
        message("Only intersecting file names were retained ...")
        message("It seems that the following genome assembly files do not have corresponding LTRpred.meta() results: ",
                paste0(dplyr::setdiff(assembly_files_chop, ltrpred_meta_results_chop), collapse = ", "))
    }
        
    files_intersect <- dplyr::intersect(assembly_files_chop, ltrpred_meta_results_chop)
    
    if (length(files_intersect) == 0)
        stop("It seems that none of the genome assembly file names in '", genome.folder, 
             "' match with the LTRpred.meta() results folder names in '",ltrpred.meta.folder,"'.",
             "Please make sure that file names match: e.g. Hsapiens.fa (assembly name in genome.folder) and Hsapiens_ltrpred (LTRpred.meta() result in ltrpred.meta.folder).",
             call. = FALSE)
    
    message("\n")
    message("Retrieving genome summaries for species: ")
    message(paste(files_intersect, collapse = ", "))
    message("\n")
    
    
    ltr_similarity_binned <- vector("list", length(files_intersect))
    n_ltrs <- vector("numeric", length(files_intersect))
    n_ltrs_freq <- vector("numeric", length(files_intersect))
    genome_size_nucl <- vector("numeric", length(files_intersect))
    genome_size_nucl_mbp <- vector("numeric", length(files_intersect))
    total_ltrs_nucl_mbp <- vector("numeric", length(files_intersect))
    total_ltrs_nucl_freq <- vector("numeric", length(files_intersect))
    NNN_freq <- vector("numeric", length(files_intersect))

    for (i in seq_len(length(files_intersect))) {
        
        ltrpred_data_sheet_path <- file.path(
            ltrpred.meta.folder,
            paste0(files_intersect[i], "_ltrpred"),
            paste0(files_intersect[i],"_LTRpred_DataSheet.tsv")
        )
        
        # determine the genome file name
        genome_file <- assembly_files[which(stringr::str_detect(assembly_files, files_intersect[i]))]
    
        message("Processing file: ", ltrpred_data_sheet_path)

        if (!file.exists(ltrpred_data_sheet_path)) {
            message("The file '", ltrpred_data_sheet_path, " could not be found ... therefore, analysis for species '",files_intersect[i],"' is omitted.")
        } else {
            pred <- read.ltrpred(ltrpred_data_sheet_path)
            
            if (quality.filter) {
                # try to reduce false positives by filtering for PBS and ORFs and rel #N's in TE <= 0.1
                pred <-
                    quality.filter(pred,
                                   sim = sim,
                                   n.orfs = n.orfs,
                                   strategy = strategy)
            }
            
            if (!quality.filter) {
                ltr_similarity <- NULL
                # keep all predicted LTR transposons including false positives
                pred <-
                    dplyr::filter(pred, ltr_similarity >= sim)
                message("No quality filter has been applied. Threshold: sim = ",
                        sim,
                        "%.")
            }
            
            # bin ltr similarities in e.g. 2% bins -> e.g. [70,72]; [72,74], etc
            binned.similarities <- cut(
                pred$ltr_similarity,
                rev(seq(100, sim,-cut.range)),
                include.lowest = TRUE,
                right = TRUE
            )
            
            # implement error handling here or a more dynamic approach
            # to handle different bin ranges in different organisms
            # sim.mass.summary <- dplyr::summarize(dplyr::group_by(pred,similarity), mass = sum(width) / 1000000)
            #
            # SimMatrix[i] <- list(sim.mass.summary)
            #
            
            ltr_similarity_binned[i] <- list(table(factor(
                binned.similarities,
                levels = levels(binned.similarities)
            )))
            
            # count the number of predicted LTR transposons
            n_ltrs[i] <- length(unique(pred$ID))
            # determine the total length of all LTR transposons in Mega base pairs
            total_ltrs_nucl_mbp[i] <- sum(pred$width) / 1000000
            # determine the genome size
            genome_size_nucl <-
                Biostrings::readDNAStringSet(file.path(genome.folder, genome_file))
            # compute genome size in Mega base pairs
            genome_size_nucl_mbp[i] <-
                sum(as.numeric(genome_size_nucl@ranges@width)) / 1000000
            
            # compute normalized LTR count: LTR count / genome size in Mbp
            n_ltrs_freq[i] <-
                as.numeric(length(unique(pred$ID)) / genome_size_nucl_mbp[i])
            # compute the proportion of LTR retrotransposons with the entire genome
            total_ltrs_nucl_freq[i] <- total_ltrs_nucl_mbp[i] / genome_size_nucl_mbp[i]
            
            # compute relative frequency of N's in genome: abs N / genome length
            NNN_freq[i] <-
                sum(as.numeric(Biostrings::vcountPattern("N", genome_size_nucl))) / sum(as.numeric(genome_size_nucl@ranges@width))
        }
    }
    
    GenomeInfo <- tibble::tibble(
        organism = files_intersect,
        total_ltrs_nucl_mbp = total_ltrs_nucl_mbp,
        total_ltrs_nucl_freq = total_ltrs_nucl_freq,
        n_ltrs = n_ltrs,
        n_ltrs_freq = n_ltrs_freq,
        genome_size_nucl_mbp = genome_size_nucl_mbp,
        NNN_freq = NNN_freq
    )
    
    SimMatrix <- do.call(rbind, ltr_similarity_binned)
    SimMatrix <-
        data.frame(organism = files_intersect, SimMatrix)
    
   
        # store results in working directory
        readr::write_delim(
            SimMatrix,
            path = paste0(file.name, "_SimilarityMatrix.csv"),
            col_names = TRUE,
            delim = ";"
        )
        
        readr::write_delim(
            GenomeInfo,
            path = paste0(file.name, "_GenomeInfo.csv"),
            col_names = TRUE,
            delim = ";"
        )
        
    message("Finished genome summary retrieval!")
}


        





