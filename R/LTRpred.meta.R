#' @title Perform Meta-Analyses with LTRpred
#' @description Run \code{\link{LTRpred}} on several genomes (sequentially) that are stored in a given folder.
#' @param genome.folder
#' @param output.name
#' @param cores
#' @param \dots 
#' @export

LTRpred.meta <- function(genome.folder,
                         output.name,
                         cores = 1,
                         ...) {
  
   if (cores > parallel::detectCores())
     stop("You sepcified more cores than are available on your machine. Please provide the correct number of cores.", call. = FALSE)
  
  
  
    # Setup cluster
    clust <- parallel::makeCluster(cores) 
    
    parallel::parLapply(clust,seq_len(length(genomes)), function (i) {
        
    })
    message("Finished meta analysis!")
  }

