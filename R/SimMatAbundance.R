#' @title Compute histogram shape similarity between species
#' @description For each pairwise species comparison the pattern similarity between
#' LTR age distributions is computed via \code{1 - cor(species1,species2)}.
#' @param sim.matrix a species age distribution matrix generated with \code{\link{LTRpred.meta}}.
#' @author Hajk-Georg Drost
#' @export
 
SimMatAbundance <- function(sim.matrix){
    res.mat <- matrix(NA_real_,nrow(sim.matrix),nrow(sim.matrix))
    nc <- vector("numeric",ncol(sim.matrix))
    nc <- ncol(sim.matrix)
    
    for (i in 1:nrow(sim.matrix)) {
        for (j in 1:nrow(sim.matrix)) {
         res.mat[i,j] <- 1 - cor(unlist(sim.matrix[i, 2:nc]),unlist(sim.matrix[j, 2:nc]))    
        }
    }
    rownames(res.mat) <- unlist(sim.matrix[ , "organism"])
    colnames(res.mat) <- unlist(sim.matrix[ , "organism"])
    
    return (res.mat)
}


