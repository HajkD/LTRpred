#' @title Plot Genome size vs. LTR transposon count for jumpers
#' @description Visualize the Genome size vs. LTR transposon count of LTR transposons
#' that have been classified as jumpers due to their LTR similarity specified
#' in the \code{ltr.similarity} arument.
#' @param genome.folder path to the folder storing all unmasked genomes for which \code{\link[LTRpred]{LTRpred}} based de novo LTR retrotransposon prediction shall be performed.
#' @param genome.matrix Genome matrix retruned by \code{\link[LTRpred]{LTRpred.meta}}.
#' @param ltr.similarity similarity threshold for defining LTR similarity. This criterion defines
#' whether LTR transposons are jumpers or not. LTR transposons fulfulling this threshold are considered
#' as being jumpers.
#' @param sci.name.len does the scientific name of the input organisms stored in \code{genome.matrix}
#' include one name (e.g. Ppaterns, = \code{sci.name.len = 1}) or two names (e.g. Physcomitrella_paterns, = \code{sci.name.len = 2}).
#' @param name.sep name separator if organism name in \code{genome.matrix} is for example \code{Ppaterns_ltrpred}. Then \code{name.sep = "_"} should be specified.
#' @param main main text of the plot.
#' @param ... additional arguments passed to \code{\link{plotSize}}.
#' @author Hajk-Georg Drost
#' @importFrom magrittr %>%
#' @export

plotSizeJumpers <- function(genome.folder, 
                              genome.matrix, 
                              ltr.similarity = 98,
                              sci.name.len   = 1,
                              name.sep       = "_",
                              main           = paste0("Genome size vs. LTR transp. count ( ",ltr.similarity,"% )"),
                                      ...){
  
  if (sci.name.len > 2)
    stop ("So far only scientific names of max. length 2 have been implemented!") 
  
  organism <- n <- nLTRs <- NULL
  
  meta.summarize.table <-
    LTRpred::meta.summarize(genome.folder, ltr.similarity)
  gm.jumpers <-
    meta.summarize.table %>% dplyr::group_by(organism) %>% dplyr::summarise(nLTRs = n())
  
  if (sci.name.len == 1) {
    gm.jumpers$organism <-
      unlist(sapply(gm.jumpers$organism, function(x) {
        unlist(stringr::str_split(x, name.sep))[1]
        
      }))
  }
  
  if (sci.name.len == 2) {
    gm.jumpers$organism <-
      unlist(sapply(gm.jumpers$organism, function(x) {
        paste0(unlist(stringr::str_split(x, name.sep))[1:2], collapse = "_")
        
      }))
  }
  
  gm.jumpers.joined <-
    dplyr::inner_join(dplyr::select(genome.matrix, -nLTRs), gm.jumpers, by = "organism")
  
  plotSize(genome.summary = gm.jumpers.joined, main = main, ...)
}