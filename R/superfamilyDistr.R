superfamilyDistr <- function(blast.file) {
  
  query_id <- subject_id <- evalue <- species <- superfamily <- NULL
  n_copia <- n_gypsy <- n_total <- repbase_blast_summarized_test_summary <- NULL
  
  repbase_blast <- readr::read_csv(blast.file, col_names = TRUE)
  
  repbase_blast_summarized <- dplyr::summarize(dplyr::group_by(repbase_blast, query_id, subject_id), evalue = min(evalue))
  species_names <- unlist(lapply(repbase_blast_summarized$query_id, function(x) unlist(stringr::str_split(x, "_"))[1]))
  superfamily_names <- unlist(lapply(repbase_blast_summarized$subject_id, function(x) unlist(stringr::str_split(x, "-"))[1]))
  name_df <- tibble::tibble(species = species_names, superfamily = superfamily_names)
  
  repbase_blast_summarized <- dplyr::bind_cols(repbase_blast_summarized, name_df)
  
  repbase_blast_summarized_summary <-
    dplyr::summarize(
      dplyr::group_by(repbase_blast_summarized, species),
      n_copia = sum(stringr::str_detect(superfamily, "Copia")),
      n_gypsy = sum(stringr::str_detect(superfamily, "Gypsy")),
      prop_copia = n_copia / (n_copia + n_gypsy),
      prop_gypsy = n_gypsy / (n_copia + n_gypsy),
      n_total = dplyr::n(),
      n_other = n_total - (n_copia + n_gypsy)
    )
  
  
  readr::write_excel_csv(repbase_blast_summarized_summary, paste0("repbase_blast_summarized_summary_", basename(blast.file) ,".csv"))
  
  return(repbase_blast_summarized_test_summary)
}
