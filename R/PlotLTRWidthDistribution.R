#' @title Plot the LTR width distribution
#' @description 
#' @param ltr.harvest.prediction
#' @param type
#' @param similarity.bin
#' @param min.sim
#' @author Hajk-Georg Drost
#' @details 
#' 
#' @examples 
#' \dontrun{
#' 
#' }
#' 
#' @export

PlotLTRWidthDistribution <- function(ltr.harvest.prediction, type = "boxplot", similarity.bin = NULL, min.sim = NULL){
    
    if (!is.element(type, c("boxplot","violin")))
        stop ("Please choose either type = 'boxplot' or type = 'violin'.")
    
    ltr_order <- width <- similarity <- ltr_similarity <- NULL
  
    if (is.null(similarity.bin) & is.null (min.sim)){
        
        #max.width <- max(ltr.harvest.prediction$ltr.retrotransposon[ , "ltr_similarity"])
        if (type == "boxplot"){
            
            res <- ggplot2::ggplot(ltr.harvest.prediction$ltr, ggplot2::aes(x = ltr_order , y = width, fill = similarity))  + ggplot2::geom_boxplot() + 
                ggplot2::facet_grid(. ~ similarity)  + ggplot2::labs(x = "\nLTR Position", y = "LTR length in bp\n") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1)) +
                ggplot2::theme(axis.text=ggplot2::element_text(size=18), axis.title=ggplot2::element_text(size=20,face="bold")) + ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange"))
            
#             res <- ggplot2::ggplot(ltr.harvest.prediction$ltr, ggplot2::aes(x = similarity , y = width, fill = similarity)) + ggplot2::geom_boxplot() + ggplot2::labs(x = "\nLTR % Similarity", y = "LTR length in bp\n") + 
#                 ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1)) +
#                 ggplot2::theme(axis.text=ggplot2::element_text(size=18), axis.title=ggplot2::element_text(size=20,face="bold")) 
            #+ ggplot2::scale_y_continuous(breaks = seq(0,ceiling(max.width),ceiling(max.width) * 0.1))
            
        }
        
        else if (type == "violin"){
            
            
            res <- ggplot2::ggplot(ltr.harvest.prediction$ltr, ggplot2::aes(x = ltr_order , y = width, fill = similarity))  + ggplot2::geom_violin() + 
                ggplot2::facet_grid(. ~ similarity)  + ggplot2::labs(x = "\nLTR Position", y = "LTR length in bp\n") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1)) +
                ggplot2::theme(axis.text=ggplot2::element_text(size=18), axis.title=ggplot2::element_text(size=20,face="bold")) + ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange"))
            
            
#             res <- ggplot2::ggplot(ltr.harvest.prediction$ltr, ggplot2::aes(x = similarity , y = width, fill = similarity)) + ggplot2::geom_violin(trim = FALSE) + ggplot2::labs(x = "\nLTR % Similarity", y = "LTR length in bp\n") + 
#                 ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1)) +
#                 ggplot2::theme(axis.text=ggplot2::element_text(size=18), axis.title=ggplot2::element_text(size=20,face="bold"))
            #+ ggplot2::scale_y_continuous(breaks = seq(0,ceiling(max.width),ceiling(max.width) * 0.1))
            
        }
    }
    
    if (!is.null(similarity.bin) & !is.null (min.sim)){
        
        ltr.harvest.prediction$ltr <- dplyr::filter(ltr.harvest.prediction$ltr, ltr_similarity >= min.sim)
        ltr.harvest.prediction$ltr <- dplyr::mutate(ltr.harvest.prediction$ltr, 
                                                                    similarity = cut(ltr_similarity,
                                                                                     rev(seq(100,min.sim,-similarity.bin)),
                                                                                     include.lowest = TRUE,
                                                                                     right          = TRUE))
        
        #max.width <- max(ltr.harvest.prediction$ltr.retrotransposon[ , "ltr_similarity"])
        if (type == "boxplot"){
            
            res <- ggplot2::ggplot(ltr.harvest.prediction$ltr, ggplot2::aes(x = ltr_order , y = width, fill = similarity))  + ggplot2::geom_boxplot() + 
                ggplot2::facet_grid(. ~ similarity)  + ggplot2::labs(x = "\nLTR Position", y = "LTR length in bp\n") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1)) +
                ggplot2::theme(axis.text=ggplot2::element_text(size=18), axis.title=ggplot2::element_text(size=20,face="bold")) + ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange"))
            
#             res <- ggplot2::ggplot(ltr.harvest.prediction$ltr,ltr_similarity, ggplot2::aes(x = similarity , y = width, fill = similarity)) + ggplot2::geom_boxplot() + ggplot2::labs(x = "\nLTR % Similarity", y = "LTR length in bp\n") + 
#                 ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1)) +
#                 ggplot2::theme(axis.text=ggplot2::element_text(size=18), axis.title=ggplot2::element_text(size=20,face="bold")) 
            #+ ggplot2::scale_y_continuous(breaks = seq(0,ceiling(max.width),ceiling(max.width) * 0.1))
            
        }
        
        else if (type == "violin"){
            
            res <- ggplot2::ggplot(ltr.harvest.prediction$ltr, ggplot2::aes(x = ltr_order , y = width, fill = similarity))  + ggplot2::geom_violin() + 
                ggplot2::facet_grid(. ~ similarity)  + ggplot2::labs(x = "\nLTR Position", y = "LTR length in bp\n") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1)) +
                ggplot2::theme(axis.text=ggplot2::element_text(size=18), axis.title=ggplot2::element_text(size=20,face="bold")) + ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange"))
            
            
#             res <- ggplot2::ggplot(ltr.harvest.prediction$ltr, ggplot2::aes(x = similarity , y = width, fill = similarity)) + ggplot2::geom_violin(trim = FALSE) + ggplot2::labs(x = "\nLTR % Similarity", y = "LTR length in bp\n") + 
#                 ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1)) +
#                 ggplot2::theme(axis.text=ggplot2::element_text(size=18), axis.title=ggplot2::element_text(size=20,face="bold"))
            #+ ggplot2::scale_y_continuous(breaks = seq(0,ceiling(max.width),ceiling(max.width) * 0.1))
            
        }
    }
    
    return(res)
}


