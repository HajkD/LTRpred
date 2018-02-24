#' @title plotSize for \code{generate.multi.quality.filter.meta} output
#' @description 
#' @param meta.summary.file
#' @param cor.method
#' @param type
#' @param label.organism
#' @param smooth.method
#' @param se
#' @param main
#' @param xlab
#' @param ylab
#' @param alpha
#' @param arrow_lab
#' @param text.size
#' @param label.size
#' @param check.overlap
#' @author Hajk-Georg Drost
#' @export

plotSizeMulti <-
  function(meta.summary.file,
           cor.method = "spearman",
           type = "total_ltrs_nucl_mbp",
           label.organism = TRUE,
           smooth.method = "lm",
           se = FALSE,
           main = "",
           xlab                = "LTR retrotransposon content in Mega [bp]",
           ylab                = "Genome size in Mega [bp]",
           alpha               = 0.3,
           arrow_lab           = FALSE,
           text.size           = 18,
           label.size          = 4,
           check.overlap       = TRUE
           ) {
    
    
  gm_files_combined <- dplyr::bind_rows(lapply(meta.summary.file, function(x) x$gm_file))  
    
  
  if (!is.element(type, c("total_ltrs_nucl_mbp", "n_ltrs_freq", "n_ltrs", "total_ltrs_nucl_freq")))
    stop("Please specify: type = 'total_ltrs_nucl_mbp' for total length of all TEs in Mbp; type = 'total_ltrs_nucl_freq' for proportion of TEs within entire genome in %; type = 'n_ltrs' for total number of TEs in genome; type = 'n_ltrs_freq' for total number of TEs in genome normalized by genome size in Mbp.", call. = FALSE)
  
  if (type == "n_ltrs") {
    n_ltrs <- NULL
    n_row_before <- nrow(gm_files_combined)
    gm_files_combined <- dplyr::filter(gm_files_combined, n_ltrs > 0)
    n_row_after <- nrow(gm_files_combined)
    
    if ((n_row_before - n_row_after) > 0)
      message(n_row_before - n_row_after, " species were removed, because they had 0 LTR retrotransposons.")
    
    message(
      "R_pearson(n_ltrs, genome_size_nucl_mbp) = ",
      cor(
        gm_files_combined$n_ltrs,
        gm_files_combined$genome_size_nucl_mbp,
        method = "pearson"
      ),
      "; p-val = ",
      cor.test(gm_files_combined$n_ltrs,
               gm_files_combined$genome_size_nucl_mbp,
               method = "pearson")
    )
    
    message(
      "R_spearman(n_ltrs, genome_size_nucl_mbp) = ",
      cor(
        gm_files_combined$n_ltrs,
        gm_files_combined$genome_size_nucl_mbp,
        method = "spearman"
      ),
      "; p-val = ",
      cor.test(gm_files_combined$n_ltrs,
               gm_files_combined$genome_size_nucl_mbp,
               method = "spearman")
    )
    
    message(
      "R_kendall(n_ltrs, genome_size_nucl_mbp) = ",
      cor(
        gm_files_combined$n_ltrs,
        gm_files_combined$genome_size_nucl_mbp,
        method = "kendall"
      ),
      "; p-val = ",
      cor.test(gm_files_combined$n_ltrs,
               gm_files_combined$genome_size_nucl_mbp,
               method = "kendall")
    )
    
    res <-
      ggplot2::ggplot(gm_files_combined,
                      ggplot2::aes(x = n_ltrs,
                                   y = genome_size_nucl_mbp,
                                   colour = kingdom)) + 
      ggplot2::facet_grid(. ~ kingdom) + ggsci::scale_color_lancet()
  }
    
  
  if (type == "n_ltrs_freq") {
    
    n_ltrs_freq <- NULL
    n_row_before <- nrow(gm_files_combined)
    gm_files_combined <- dplyr::filter(gm_files_combined, n_ltrs_freq > 0)
    n_row_after <- nrow(gm_files_combined)
    
    if ((n_row_before - n_row_after) > 0)
      message(n_row_before - n_row_after, " species were removed, because they had 0 LTR retrotransposons.")
    
    
    message(
      "R_pearson(n_ltrs_freq, genome_size_nucl_mbp) = ",
      cor(
        gm_files_combined$n_ltrs_freq,
        gm_files_combined$genome_size_nucl_mbp,
        method = "pearson"
      ),
      "; p-val = ",
      cor.test(gm_files_combined$n_ltrs_freq,
               gm_files_combined$genome_size_nucl_mbp,
               method = "pearson")
    )
    
    
    message(
      "R_spearman(n_ltrs_freq, genome_size_nucl_mbp) = ",
      cor(
        gm_files_combined$n_ltrs_freq,
        gm_files_combined$genome_size_nucl_mbp,
        method = "spearman"
      ),
      "; p-val = ",
      cor.test(gm_files_combined$n_ltrs_freq,
               gm_files_combined$genome_size_nucl_mbp,
               method = "spearman")
    )
    
    message(
      "R_kendall(n_ltrs_freq, genome_size_nucl_mbp) = ",
      cor(
        gm_files_combined$n_ltrs_freq,
        gm_files_combined$genome_size_nucl_mbp,
        method = "kendall"
      ),
      "; p-val = ",
      cor.test(gm_files_combined$n_ltrs_freq,
               gm_files_combined$genome_size_nucl_mbp,
               method = "kendall")
    )
    
    res <-
      ggplot2::ggplot(gm_files_combined,
                      ggplot2::aes(x = n_ltrs_freq,
                                   y = genome_size_nucl_mbp,
                                   colour = kingdom)) + 
      ggplot2::facet_grid(. ~ kingdom) + ggsci::scale_color_lancet()
  }
    
  
  if (type == "total_ltrs_nucl_mbp") {
    total_ltrs_nucl_mbp <- NULL
    n_row_before <- nrow(gm_files_combined)
    gm_files_combined <- dplyr::filter(gm_files_combined, total_ltrs_nucl_mbp > 0)
    n_row_after <- nrow(gm_files_combined)
    
    if ((n_row_before - n_row_after) > 0)
      message(n_row_before - n_row_after, " species were removed, because they had 0 LTR retrotransposons.")
    
    print(dplyr::summarise(dplyr::group_by(gm_files_combined, kingdom), 
      R_pearson = round(cor(
      total_ltrs_nucl_mbp,
      genome_size_nucl_mbp,
      method = "pearson"
    ), 3),
    p_val_pearson = cor.test(
      total_ltrs_nucl_mbp,
      genome_size_nucl_mbp,
      method = "pearson"
    )$p.value,
    R_spearman = round(cor(
      total_ltrs_nucl_mbp,
      genome_size_nucl_mbp,
      method = "spearman"
    ), 3),
    p_val_spearman = cor.test(
      total_ltrs_nucl_mbp,
      genome_size_nucl_mbp,
      method = "spearman"
    )$p.value,
    R_kendall = round(cor(
      total_ltrs_nucl_mbp,
      genome_size_nucl_mbp,
      method = "kendall"
    ), 3),
    p_val_kendall = cor.test(
      total_ltrs_nucl_mbp,
      genome_size_nucl_mbp,
      method = "kendall"
    )$p.value
    ))
    
    res <-
      ggplot2::ggplot(
        gm_files_combined,
        ggplot2::aes(x = total_ltrs_nucl_mbp,
                     y = genome_size_nucl_mbp,
                     colour = kingdom)
      ) + ggplot2::facet_grid(. ~ kingdom) + ggsci::scale_color_lancet()
  }
    
  
  if (type == "total_ltrs_nucl_freq") {
    total_ltrs_nucl_freq <- NULL
    n_row_before <- nrow(gm_files_combined)
    gm_files_combined <- dplyr::filter(gm_files_combined, total_ltrs_nucl_freq > 0)
    n_row_after <- nrow(gm_files_combined)
    
    if ((n_row_before - n_row_after) > 0)
      message(n_row_before - n_row_after, " species were removed, because they had 0 LTR retrotransposons.")
    
    message(
      "R_spearman(total_ltrs_nucl_freq, genome_size_nucl_mbp) = ",
      cor(
        gm_files_combined$total_ltrs_nucl_freq,
        gm_files_combined$genome_size_nucl_mbp,
        method = "spearman"
      ),
      "; p-val = ",
      cor.test(m_files_combined$total_ltrs_nucl_freq,
               gm_files_combined$genome_size_nucl_mbp,
               method = "spearman")
    )
    
    message(
      "R_kendall(total_ltrs_nucl_freq, genome_size_nucl_mbp) = ",
      cor(
        gm_files_combined$total_ltrs_nucl_freq,
        gm_files_combined$genome_size_nucl_mbp,
        method = "kendall"
      ),
      "; p-val = ",
      cor.test(m_files_combined$total_ltrs_nucl_freq,
               gm_files_combined$genome_size_nucl_mbp,
               method = "kendall")
    )
    
    
    res <-
      ggplot2::ggplot(
        gm_files_combined,
        ggplot2::aes(x = total_ltrs_nucl_freq * 100,
                     y = genome_size_nucl_mbp,
                     colour = kingdom)
      ) + ggplot2::facet_grid(. ~ kingdom) + 
      ggsci::scale_color_lancet()
  }
    
  
  res <- res + ggplot2::geom_point(size = 3, show.legend = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = xlab,
      y = ylab) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = text.size)) +
    ggplot2::theme(
      axis.title  = ggplot2::element_text(size = text.size, face = "bold"),
      axis.text.y = ggplot2::element_text(size = text.size, face = "bold"),
      axis.text.x = ggplot2::element_text(size = text.size, face = "bold"),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        size = text.size,
        colour = "black",
        face = "bold"
      )
    )
  
  if (label.organism) {
    if (!arrow_lab) {
      res <- res + ggplot2::geom_text(
        ggplot2::aes(label = organism),
        hjust = 0,
        vjust = 0,
        size = label.size,
        check_overlap = check.overlap,
        fontface = 'bold'
      )
    }
    
    if (arrow_lab) {
      res <-
        res + ggrepel::geom_text_repel(ggplot2::aes(label = organism),
                                       size = label.size,
                                       fontface = 'bold')
    }
  }
  
  res <- res +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6))
  
  if (!is.null(smooth.method)) {
    res <-
      res + ggplot2::geom_smooth(method = smooth.method,
                                 se = se,
                                 alpha = alpha,
                                 colour = "black",
                                 size = 2)
  }
  
  return(res)

}
