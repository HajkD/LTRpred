#' @title Plot Ranges of an genomic feature
#' @description Helper function to plot the genomic range of the predicted LTR transposon.
#' @param x range object.
#' @param xlim range of the x-axis.
#' @param main main text for the plot generated. 
#' @param col color of the visualized ranges.
#' @param sep separation between ranges.
#' @param ... standard graphics parameters passed to \code{link{rect}}.
#' @author Hajk-Georg Drost


PlotRanges <- function(x, 
                       xlim = x, 
                       main = deparse(substitute(x)), 
                       col  = "black", 
                       sep  = 0.5, ...){
    height <- 1
    
    if (is(xlim, "Ranges"))
        xlim <- c(0, max(IRanges::end(xlim)))
    
    bins <- IRanges::disjointBins(IRanges::IRanges(IRanges::start(x), IRanges::end(x) + 1))
    
    graphics::plot.new()
    
    graphics::plot.window(xlim, c(0, max(bins)*(height + sep)))
    
    ybottom <- bins * (sep + height) - height
    
    graphics::rect(IRanges::start(x)-0.5, ybottom, IRanges::end(x) + 0.5, ybottom + height, col = col, ...)
    
    graphics::title(main)
    graphics::axis(1) 
}



PlotRanges2 <- function(x, xlim = x, main = deparse(substitute(x)), col = "black", sep = 0.5, ...){
    height <- 1
    
    if (methods::is(xlim, "Ranges"))
        xlim <- c(min(IRanges::start(xlim)), max(IRanges::end(xlim)))
    
    bins <- IRanges::disjointBins(IRanges::IRanges(IRanges::start(x), IRanges::end(x) + 1))
    
    graphics::plot.new()
    
    graphics::plot.window(xlim, c(0, max(bins)*(height + sep)))
    
    ybottom <- bins * (sep + height) - height
    
    graphics::rect(IRanges::start(x)-0.5, ybottom, IRanges::end(x) + 0.5, ybottom + height, col = col, ...)
    
    graphics::title(main)
    graphics::axis(1) 
}











