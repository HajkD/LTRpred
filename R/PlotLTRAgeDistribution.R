#' @title Plot the age distribution of predicted LTR transposons
#' @description This function
#' @param ltr.harvest.prediction
#' @param similarity.bin
#' @param min.sim
#' @author Hajk-Georg Drost
#' @details This
#' 
#' @examples 
#' \dontrun{
#' 
#' }
#' 
#' @export

PlotLTRAgeDistribution <- function(ltr.harvest.prediction, similarity.bin = NULL, min.sim = NULL){
    
    similarity <- ltr_similarity <- NULL
    
    if (is.null(similarity.bin) & is.null (min.sim)){
        
        res <- ggplot2::ggplot(ltr.harvest.prediction, ggplot2::aes(x = similarity), order = FALSE) +
            ggplot2::geom_bar(stat="bin", fill = "turquoise4") + ggplot2::labs(x = "LTR % Similarity", y = "Frequency") +
            ggplot2::theme_minimal() + ggplot2::theme(axis.title=ggplot2::element_text(size=14,face="bold"),axis.text.y = ggplot2::element_text(size=14,face="bold"),axis.text.x = ggplot2::element_text(size=14,face="bold"),
                                                      panel.background = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 20, colour = "black", face = "bold"))
        
    }
    
    if (!is.null(similarity.bin) & !is.null (min.sim)){
        
        ltr.harvest.prediction <- dplyr::filter(ltr.harvest.prediction, ltr_similarity >= min.sim)
        ltr.harvest.prediction <- dplyr::mutate(ltr.harvest.prediction, 
                                                    similarity = cut(ltr_similarity,
                                                                     rev(seq(100,min.sim,-similarity.bin)),
                                                                     include.lowest = TRUE,
                                                                     right          = TRUE))
        
    
        res <- ggplot2::ggplot(ltr.harvest.prediction, ggplot2::aes(x = similarity), order = FALSE) +
        ggplot2::geom_bar(stat="bin", fill = "turquoise4") + ggplot2::labs(x = "LTR % Similarity", y = "Frequency") +
        ggplot2::theme_minimal() + ggplot2::theme(axis.title=ggplot2::element_text(size=14,face="bold"),axis.text.y = ggplot2::element_text(size=14,face="bold"),axis.text.x = ggplot2::element_text(size=14,face="bold"),
        panel.background = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 20, colour = "black", face = "bold"))
    }
    
    return(res)
}