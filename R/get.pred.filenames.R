#' @title Retrieve file names of files genereated by LTRpred
#' @description Useful helper function to generate automated file names
#' when crawling LTRpred generated result folders.
#' @param LTRpred.folder path to a \code{\link{LTRpred}} generated folder.
#' @param ws.wrap shall white space separated path names be wrapped with \code{'file name'}.
#' @author Hajk-Georg Drost
#' @examples 
#' # without wrapped whitespaces
#' get.pred.filenames("path/to/LTRpred.folder")
#' 
#' # with wrapped whitespaces
#' get.pred.filenames("path to/LTRpred.folder", ws.wrap = TRUE )
#' @export

get.pred.filenames <- function(LTRpred.folder, ws.wrap = FALSE){
  
  folder.name <- basename(LTRpred.folder)
  folder.name <- stringr::str_replace(folder.name,"_ltrpred","")
  
  ltrdigest.folder <- paste0(folder.name,"_ltrdigest")
  
  pred.files <- c(
    file.path(LTRpred.folder,paste0(folder.name,"_LTRpred_DataSheet.csv")),
    file.path(LTRpred.folder,paste0(folder.name,"_LTRpred.bed")),
    file.path(LTRpred.folder,paste0(folder.name,"_LTRpred.gff")),
    file.path(LTRpred.folder,paste0(folder.name,"_orfs_nt.fsa")),
    file.path(LTRpred.folder,paste0(folder.name,"_orfs_aa.fsa")),
    file.path(LTRpred.folder,paste0(folder.name,"-ltrdigest_complete.fas_nt.fsa")),
    file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_3ltr.fas")),
    file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_5ltr.fas")),
    file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_5ltr.fas")),
    file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"_LTRdigestPrediction.gff")),
    file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_tabout.csv")),
    file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_complete.fas")),
    file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_conditions.csv")),
    file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_pbs.fas")),
    file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_ppt.fas")),
    file.path(LTRpred.folder,ltrdigest.folder,paste0(folder.name,"-ltrdigest_pbs.fas")))
  
  if (ws.wrap){
    pred.files <- ws.wrap.path(pred.files)
  }
  
 return(pred.files)
}