#' @title Meta function for applying \code{rename.fasta}
#' @description Meta function for applying \code{\link{rename.fasta}} to many fasta files within a folder.
#' @param in.folder path to folder storing \code{\link{LTRpred}} generated organism prediction folders.
#' @param out.file name of concatenated output file. 
#' @author Hajk-Georg Drost
#' @export

rename.fasta.meta <- function(in.folder, out.file){
    
    if (!file.exists(in.folder))
        stop("The folder '",
             in.folder,
             "' does not exist. Please provide a correct path.")
    
    LTRpred.folders <- list.files(in.folder)
    
    for (i in seq_len(length(LTRpred.folders))) {
      message("Processing file ",LTRpred.folders[i]," ...")
        species.name <-
            stringr::str_replace(LTRpred.folders[i], "_ltrpred", "")
        rename.fasta(
            file.path(
                in.folder,
                LTRpred.folders[i],
                paste0(species.name, "_ltrdigest"),
                paste0(species.name, "-ltrdigest_complete.fas")
            ),
            species = species.name,
            output = out.file,
            append = TRUE
        )
    }
}
