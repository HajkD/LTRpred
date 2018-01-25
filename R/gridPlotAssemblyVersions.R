#' @title Plot retrotransposon age distributions for predictions of different genome assembly versions
#' @description Compare the abundance of retrotransposon predictions between different genome assembly versions. This way a genome assembly version bias can be either ruled out
#' or considered for subsequent conclusions.
#' @param data prediction files returned by \code{\link[LTRpred]{LTRpred}}. Usually a combined prediction file that was generated with \code{\link{combinePreds}}.
#' @param quality.filter shall false positives be filtered out as much as possible or not.
#' @param text.size size of x-axis, y-axis, and title text.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param main.text title text.
#' @param y.ticks number of ticks that shall be drawn on the y-axis.
#' @param sim If \code{quality.filter = TRUE}: LTR similarity threshold. Only putative LTR transposons that fulfill this 
#' LTR similarity threshold will be retained.
#' @param n.orfs If \code{quality.filter = TRUE}: minimum number of ORFs detected in the putative LTR transposon.
#' @author Hajk-Georg Drost
#' @export
gridPlotAssemblyVersions <-
  function(data,
           quality.filter = FALSE,
           text.size = 18,
           xlab = "LTR similarity",
           ylab = "LTR retrotransposon count",
           main.text = "",
           y.ticks = 6,
           sim = 70,
           n.orfs = 0) {
    
    similarity <- NULL
    p <-
      ggplot2::ggplot(data, ggplot2::aes(
        x = factor(similarity, levels = unique(as.character(
          cut(
            seq(70, 100, 1),
            rev(seq(100, 70,-2)),
            include.lowest = TRUE,
            right          = TRUE
          )
        ))),
        fill = factor(similarity, levels = unique(as.character(
          cut(
            seq(70, 100, 1),
            rev(seq(100, 70,-2)),
            include.lowest = TRUE,
            right          = TRUE
          )
        )))
      )) +
      ggplot2::geom_bar(stat = "count")  +
      ggplot2::facet_grid(. ~ factor(species)) +
      ggplot2::scale_fill_discrete(name = "similarity")  +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = xlab,
                    y = ylab,
                    title = main.text) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = text.size)) +
      ggplot2::theme(
        axis.title  = ggplot2::element_text(size = text.size, face = "bold"),
        axis.text.y = ggplot2::element_text(size = text.size, face = "bold"),
        axis.text.x = ggplot2::element_text(size = text.size, face = "bold"),
        strip.text.x = ggplot2::element_text(size = 14, face = "bold"),
        panel.background = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          size = text.size,
          colour = "black",
          face = "bold"
        )
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 1,
        hjust = 1
      )) +
      ggplot2::scale_fill_manual(values = customColors(length(unique(
        cut(
          seq(70, 100, 1),
          rev(seq(100, 70,-2)),
          include.lowest = TRUE,
          right          = TRUE
        )
      ))), name = "Similarity") +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y.ticks))
    return(p)
  }