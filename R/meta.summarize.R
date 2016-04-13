meta.summarize <- function(result.folder, ltr.similarity = 98){
  
  result.files <- list.files(result.folder)
  folders0 <- result.files[stringr::str_detect(result.files, "ltrpred")]
  org.list <- vector("list",length(folders0))
  ltr_similarity <- NULL
  
  for (i in 1:length(folders0)){
    choppedFolder <- unlist(stringr::str_split(folders0[i],"_"))
    pred <- readr::read_delim(file.path(result.folder,folders0[i],paste0(paste0(choppedFolder[-length(choppedFolder)],collapse = "_"),"_LTRpred_DataSheet.csv")), delim = ";")
    
    pred.filtered <- dplyr::filter(pred, ltr_similarity >= ltr.similarity)
    org.list[i] <- list(data.frame(organism = rep(stringr::str_replace(folders0[i],"_ltrpred",""),nrow(pred.filtered)),as.data.frame(pred.filtered)))
  }
  
  return (do.call(rbind,org.list))
}