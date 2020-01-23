#' @title plotSize for \code{generate.multi.quality.filter.meta} output
#' @description Visualize the correlation between genome size and TE load for \code{\link{generate.multi.quality.filter.meta}} output.
#' @param meta.summary.file a meta.summary.file file.
#' @param cor.method correlation method.
#' @param type type of TE load metric that shall be visualized. Options are
#' \itemize{
#' \item \code{type = "total_ltrs_nucl_mbp"}
#' \item \code{type = "n_ltrs_freq"}
#' \item \code{type = "n_ltrs"}
#' \item \code{type = "total_ltrs_nucl_freq"}
#' }
#' @param label.organism shall organism names be labeled next to each dot in the correlation plot?
#' @param smooth.method shall a smoothing function be applied to the correlation plot?
#' @param se shall standard.error be drawn when \code{smooth.method = TRUE}.
#' @param main title text.
#' @param xlab x-axis text.
#' @param ylab y-axis text.
#' @param alpha alpha of standard error band.
#' @param arrow_lab shall arrows be drawn between dots and label names?
#' @param text.size general text size.
#' @param label.size label text size.
#' @param check.overlap apply overlap correction for label names.
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
           xlab = "LTR retrotransposon content in Mega [bp]",
           ylab = "Genome size in Mega [bp]",
           alpha = 0.3,
           arrow_lab = FALSE,
           text.size = 18,
           label.size = 4,
           check.overlap = TRUE
           ) {
    
    n_ltrs <- total_ltrs_nucl_mbp <- 
    
  gm_files_combined <- dplyr::bind_rows(lapply(meta.summary.file, function(x) x$gm_file))  
    
  if (!is.element(type, c("total_ltrs_nucl_mbp", "n_ltrs_freq", "n_ltrs", "total_ltrs_nucl_freq")))
    stop("Please specify: type = 'total_ltrs_nucl_mbp' for total length of all TEs in Mbp; type = 'total_ltrs_nucl_freq' for proportion of TEs within entire genome in %; type = 'n_ltrs' for total number of TEs in genome; type = 'n_ltrs_freq' for total number of TEs in genome normalized by genome size in Mbp.", call. = FALSE)
  
  if (type == "n_ltrs") {
    
    
    message(
      "R_pearson(n_ltrs, genome_size_nucl_mbp) = ",
      stats::cor(
        gm_files_combined$n_ltrs,
        gm_files_combined$genome_size_nucl_mbp,
        method = "pearson"
      ),
      "; p-val = ",
      stats::cor.test(gm_files_combined$n_ltrs,
               gm_files_combined$genome_size_nucl_mbp,
               method = "pearson")
    )
    
    message(
      "R_spearman(n_ltrs, genome_size_nucl_mbp) = ",
      stats::cor(
        gm_files_combined$n_ltrs,
        gm_files_combined$genome_size_nucl_mbp,
        method = "spearman"
      ),
      "; p-val = ",
      stats::cor.test(gm_files_combined$n_ltrs,
               gm_files_combined$genome_size_nucl_mbp,
               method = "spearman")
    )
    
    message(
      "R_kendall(n_ltrs, genome_size_nucl_mbp) = ",
      stats::cor(
        gm_files_combined$n_ltrs,
        gm_files_combined$genome_size_nucl_mbp,
        method = "kendall"
      ),
      "; p-val = ",
      stats::cor.test(gm_files_combined$n_ltrs,
               gm_files_combined$genome_size_nucl_mbp,
               method = "kendall")
    )
    
    res <-
      ggplot2::ggplot(gm_files_combined,
                      ggplot2::aes(x = n_ltrs,
                                   y = genome_size_nucl_mbp,
                                   colour = kingdom)) + 
      ggplot2::facet_grid(. ~ kingdom, scales = "free_x") + ggsci::scale_color_lancet()
  }
    
  
  if (type == "n_ltrs_freq") {
  
    
    message(
      "R_pearson(n_ltrs_freq, genome_size_nucl_mbp) = ",
      stats::cor(
        gm_files_combined$n_ltrs_freq,
        gm_files_combined$genome_size_nucl_mbp,
        method = "pearson"
      ),
      "; p-val = ",
      stats::cor.test(gm_files_combined$n_ltrs_freq,
               gm_files_combined$genome_size_nucl_mbp,
               method = "pearson")
    )
    
    message(
      "R_spearman(n_ltrs_freq, genome_size_nucl_mbp) = ",
      stats::cor(
        gm_files_combined$n_ltrs_freq,
        gm_files_combined$genome_size_nucl_mbp,
        method = "spearman"
      ),
      "; p-val = ",
      stats::cor.test(gm_files_combined$n_ltrs_freq,
               gm_files_combined$genome_size_nucl_mbp,
               method = "spearman")
    )
    
    message(
      "R_kendall(n_ltrs_freq, genome_size_nucl_mbp) = ",
      stats::cor(
        gm_files_combined$n_ltrs_freq,
        gm_files_combined$genome_size_nucl_mbp,
        method = "kendall"
      ),
      "; p-val = ",
      stats::cor.test(gm_files_combined$n_ltrs_freq,
               gm_files_combined$genome_size_nucl_mbp,
               method = "kendall")
    )
    
    
    n_ltrs_freq <- genome_size_nucl_mbp <- kingdom <- NULL
      
    res <-
      ggplot2::ggplot(gm_files_combined,
                      ggplot2::aes(x = n_ltrs_freq,
                                   y = genome_size_nucl_mbp,
                                   colour = kingdom)) + 
      ggplot2::facet_grid(. ~ kingdom) + ggsci::scale_color_lancet()
  }
    
  
  if (type == "total_ltrs_nucl_mbp") {
    
    cor_data <- dplyr::summarise(dplyr::group_by(gm_files_combined, kingdom), 
      R_pearson = round(stats::cor(
      total_ltrs_nucl_mbp,
      genome_size_nucl_mbp,
      method = "pearson"
    ), 3),
    p_val_pearson = stats::cor.test(
      total_ltrs_nucl_mbp,
      genome_size_nucl_mbp,
      method = "pearson"
    )$p.value,
    R_spearman = round(stats::cor(
      total_ltrs_nucl_mbp,
      genome_size_nucl_mbp,
      method = "spearman"
    ), 3),
    p_val_spearman = stats::cor.test(
      total_ltrs_nucl_mbp,
      genome_size_nucl_mbp,
      method = "spearman"
    )$p.value,
    R_kendall = round(stats::cor(
      total_ltrs_nucl_mbp,
      genome_size_nucl_mbp,
      method = "kendall"
    ), 3),
    p_val_kendall = stats::cor.test(
      total_ltrs_nucl_mbp,
      genome_size_nucl_mbp,
      method = "kendall"
    )$p.value
    )
    
    print(cor_data)
    
    kingdom_name <-
      unique(unlist(lapply(unique(gm_files_combined$kingdom), function(x)
        unlist(stringr::str_split(x, "_"))[1])))
    
    readr::write_excel_csv(cor_data,
                           paste0("cor_data_", kingdom_name , ".csv"),
                           col_names = TRUE)
    
    res <-
      ggplot2::ggplot(
        gm_files_combined,
        ggplot2::aes(x = total_ltrs_nucl_mbp,
                     y = genome_size_nucl_mbp,
                     colour = kingdom)
      ) + ggplot2::facet_grid(. ~ kingdom, scales = "free_x") + ggsci::scale_color_lancet()
  }
    
  
  if (type == "total_ltrs_nucl_freq") {
    
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
      stats::cor(
        gm_files_combined$total_ltrs_nucl_freq,
        gm_files_combined$genome_size_nucl_mbp,
        method = "kendall"
      ),
      "; p-val = ",
      stats::cor.test(m_files_combined$total_ltrs_nucl_freq,
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
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8))
  
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
