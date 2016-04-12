#' @title Perform Meta-Analyses with LTRpred
#' @description Run \code{\link{LTRpred}} on several genomes (sequentially) that are stored in a given folder.
#' @param genome.folder path to the folder storing all unmasked genomes for which \code{\link{LTRpred}}
#' based \code{de novo} LTR retrotransposon prediction shall be performed.
#' @param result.folder folder in which \code{LTRpred} results shall be stored.
#' @param similarity similarity threshold for defining LTR similarity.
#' @param LTRpred.meta.folder meta-folder storing already pre-cumputed \code{\link{LTRpred}} generated files.
#' @param \dots all parameters of \code{\link{LTRpred}}.
#' @author Hajk-Georg Drost
#' @details 
#' This function provides a crawler to run \code{\link{LTRpred}} sequencially
#' on any number of genomes stored within the same folder.
#' The result will be saved on in the \code{result.folder} and can be analysed with the
#' \code{\link{meta.apply}} function.
#' @examples 
#' \dontrun{
#' # perform a meta analysis on multiple genomes
#' # stored in a genomes folder
#' LTRpred.meta(genome.folder = "Genomes/",
#'              result.folder = "LTRpredResults"
#'              trnas         = "plantRNA_Arabidopsis.fsa",
#'              hmms          = "hmm_*")
#' }
#' @export
   
LTRpred.meta <- function(genome.folder       = NULL, 
                         result.folder       = NULL,
                         similarity          = 70,
                         LTRpred.meta.folder = NULL, 
                         ...){
  
  
  if (!is.null(genome.folder) && is.null(result.folder) && !is.null(LTRpred.meta.folder)){
    
    if (!file.exists(genome.folder))
      stop ("The folder ' ",genome.folder," ' could not be found.")
    
    cat("\n")
    cat("Starting LTRpred meta analysis on the following files: ")
    genomes <- list.files(genome.folder)
    cat("\n")
    cat("\n")
    cat(paste(list.files(LTRpred.meta.folder), collapse = ", "))
    cat("\n")
    cat("\n")
    
    result.files <- list.files(LTRpred.meta.folder)
    folders0 <- result.files[stringr::str_detect(result.files, "ltrpred")]
    
    genomes.chopped <- sapply(genomes, function(x) unlist(stringr::str_split(x,"[.]"))[1])
    
    ltrdigest.folder.files.chopped <- sapply(folders0, function(x) unlist(stringr::str_replace(x,"_ltrpred","")))
    
    available.genomes <- match(ltrdigest.folder.files.chopped,genomes.chopped)
    genomes <- genomes[available.genomes]
    
    if (length(folders0) != length(genomes))
      stop ("Please make sure that the number of your genome files matches with your LTRpred folders.")
    
    SimMatrix <- 1
    nLTRs <- 1
    nLTRs.normalized <- 1
    gs <- 1
    
    cat("Processing file:")
    cat("\n")
    for (i in 1:length(folders0)){
      choppedFolder <- unlist(stringr::str_split(folders0[i],"_"))
      print(file.path(LTRpred.meta.folder,folders0[i],paste0(paste0(choppedFolder[-length(choppedFolder)],collapse = "_"),"_LTRpred_DataSheet.csv")))
      cat("\n")
      cat("\n")
      
      if (!file.exists(file.path(LTRpred.meta.folder,folders0[i],paste0(paste0(choppedFolder[-length(choppedFolder)],collapse = "_"),"_LTRpred_DataSheet.csv")))){
        
        print(paste0("Skip :",file.path(LTRpred.meta.folder,folders0[i],paste0(paste0(choppedFolder[-length(choppedFolder)],collapse = "_"),"_LTRpred_DataSheet.csv")), " -> folder was empty!"))
        cat("\n")
      } else {
        
        pred <- readr::read_delim(file.path(LTRpred.meta.folder,folders0[i],paste0(paste0(choppedFolder[-length(choppedFolder)],collapse = "_"),"_LTRpred_DataSheet.csv")), delim = ";")
        
        SimMatrix[i] <- list(table(factor(pred$similarity, levels = levels(cut(pred$ltr_similarity, rev(seq(100,similarity,-2)),include.lowest = TRUE,right = TRUE)))))
        nLTRs[i] <- nrow(pred)
        genome.size <-  Biostrings::readDNAStringSet(file.path(genome.folder,genomes[i]))
        gs[i] <- sum(genome.size@ranges@width)
        nLTRs.normalized[i] <- length(unique(pred$ID)) / sum(genome.size@ranges@width)
      }
}
      
    names(nLTRs) <- folders0
    names(nLTRs.normalized) <- folders0
    names(gs) <- folders0
    
    GenomeInfo <- data.frame(organism = genomes, nLTRs = nLTRs, norm.nLTRs = nLTRs.normalized, genome.size = gs )
    SimMatrix <- do.call(rbind, SimMatrix)
    SimMatrix <- data.frame(organism = folders0 , SimMatrix)
    write.table(SimMatrix,paste0(basename(LTRpred.meta.folder),"_SimilarityMatrix.csv"), sep = ";", quote = FALSE, col.names = TRUE, row.names = FALSE)
    write.table(GenomeInfo,paste0(basename(LTRpred.meta.folder),"_GenomeInfo.csv"), sep = ";", quote = FALSE, col.names = TRUE, row.names = FALSE)
    
    cat("Finished meta analysis!")
    
  } else {
    
    if (!file.exists(genome.folder))
      stop ("The folder ' ",genome.folder," ' could not be found.")
    
    cat("\n")
    cat("Starting LTRpred meta analysis on the following genomes: ")
    genomes <- list.files(genome.folder)
    cat("\n")
    cat("\n")
    cat(paste(genomes, collapse = ", "))
    cat("\n")
    cat("\n")
    genome.names.chopped <- sapply(genomes, function(x) unlist(stringr::str_split(x, "[.]"))[1])
    
    # run meta analysis for all species sequencially
    for (i in 1:length(genomes)){
      LTRpred(genome.file = file.path(genome.folder, genomes[i]), ...)
    }
    
    if (!is.null(result.folder)){
      # store results in result folder -> default: working directory
      file.move(from = paste0(genome.names.chopped,"_ltrpred"), 
                to   = file.path(result.folder,paste0(genome.names.chopped,"_ltrpred")))
    }
    
    result.files <- list.files(result.folder)
    folders0 <- result.files[stringr::str_detect(result.files, "ltrpred")]
    SimMatrix <- 1
    nLTRs <- 1
    nLTRs.normalized <- 1
    gs <- 1
    for (i in 1:length(folders0)){
      choppedFolder <- unlist(stringr::str_split(folders0[i],"_"))
      pred <- readr::read_delim(file.path(result.folder,folders0[i],paste0(paste0(choppedFolder[-length(choppedFolder)],collapse = "_"),"_LTRpred_DataSheet.csv")), delim = ";")
      
      SimMatrix[i] <- list(table(factor(pred$similarity, levels = levels(cut(pred$ltr_similarity, rev(seq(100,similarity,-2)),include.lowest = TRUE,right = TRUE)))))
      nLTRs[i] <- nrow(pred)
      genome.size <-  Biostrings::readDNAStringSet(file.path(genome.folder,genomes[i]))
      gs[i] <- sum(genome.size@ranges@width)
      nLTRs.normalized[i] <- length(unique(pred$ID)) / sum(genome.size@ranges@width)
    }
    
    names(nLTRs) <- folders0
    names(nLTRs.normalized) <- folders0
    names(gs) <- folders0
    
    GenomeInfo <- data.frame(organism = genomes, nLTRs = nLTRs, norm.nLTRs = nLTRs.normalized, genome.size = gs )
    SimMatrix <- do.call(rbind, SimMatrix)
    SimMatrix <- data.frame(organism = folders0 , SimMatrix)
    write.table(SimMatrix,paste0(basename(result.folder),"_SimilarityMatrix.csv"), sep = ";", quote = FALSE, col.names = TRUE, row.names = FALSE)
    write.table(GenomeInfo,paste0(basename(result.folder),"_GenomeInfo.csv"), sep = ";", quote = FALSE, col.names = TRUE, row.names = FALSE)
    
    cat("Finished meta analysis!")
  }
}








