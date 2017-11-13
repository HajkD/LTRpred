#' @title Colors scheme for plots
#' @description A distinguished color scheme for scientific plots.
#' @param n number of different colors that shall be returned.
#' @author Hajk-Georg Drost
#' @examples
#' # return 5 colors 
#' bcolor(5)
#' @export
bcolor <- function(n) return(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(n))
