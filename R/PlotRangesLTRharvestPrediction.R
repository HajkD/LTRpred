#' @title Plot the location of putative LTR transposons along the chromosomes
#' @description This function 
#' @param ltr.harvest.prediction
#' @author Hajk-Georg Drost
#' @export

PlotRangesLTRharvestPrediction <- function(ltr.harvest.prediction){
    
    chromosome <- NULL
    chromosomes <- names(table(ltr.harvest.prediction[ , "chromosome"]))
    par(mfrow = n2mfrow(length(chromosomes)))
    sapply(chromosomes, function(chr){
        
        ltr.ranges <- dplyr::filter(ltr.harvest.prediction, chromosome == chr)
        PlotRanges(IRanges::IRanges(ltr.ranges[ , "start"],ltr.ranges[ , "end"]), main = chr)
    })
    
}

PlotRangesLTRharvestPrediction2 <- function(ltr.harvest.prediction){
    
    chromosome <- NULL
    chromosomes <- names(table(ltr.harvest.prediction$ltr.retrotransposon[ , "chromosome"]))
#     par(mfrow = n2mfrow(length(chromosomes)))
#     sapply(chromosomes, function(chr){
#         
#         ltr.ranges <- dplyr::filter(ltr.harvest.prediction$ltr.retrotransposon, chromosome == chr)
#         ggbio::autoplot(IRanges::IRanges(ltr.ranges[ , "start"],ltr.ranges[ , "end"]), main = chr, geom = "segment")
#     })
#     
    
    ltr.ranges <- dplyr::filter(ltr.harvest.prediction$ltr.retrotransposon, chromosome == chromosomes[1])
    ggbio::autoplot(IRanges::IRanges(ltr.ranges[ , "start"],ltr.ranges[ , "end"]), main = chromosomes[1], geom = "arrowrect")
}
