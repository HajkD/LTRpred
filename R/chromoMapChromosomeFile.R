#' @title Generate chromsome file for chromoMap package visualization
#' @description Generate chromsome file for chromoMap package visualization.
#' @param data a \code{\link{data.frame}} with the following column format:
#'  \itemize{
#'  \item column 1: \code{data$chromosome_name}
#'  \item column 2: \code{data$chromosome_start}
#'  \item column 3: \code{data$chromosome_end}
#'  \item column 4: \code{data$centromere_start}
#'  }
#' @param output a file path (file name) to the output file.
#' @author Hajk-Georg Drost
#' @references https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html
#' @export   
chromoMapChromosomeFile <- function(data, output) {
  
  res <- dplyr::data_frame(chromosome_name = data$chromosome_name,
                           chromosome_start = data$chromosome_start, 
                           chromosome_end      = data$chromosome_end, 
                           centromere_start = data$centromere_start)
  
  utils::write.table(res, output, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
}