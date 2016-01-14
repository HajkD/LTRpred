#' @title Import the output of LTRharvest or LTRdigest
#' @description This function
#' @param gff.file
#' @param tabout.file
#' @param program
#' @param ltr.fasta
#' @param inner.seq.fasta
#' @param data
#' @param similarity.bin
#' @param min.sim
#' @author Hajk-Georg Drost
#' @examples 
#' \dontrun{
#' 
#' }
#' @export

read.prediction <- function( gff.file        = NULL,
                             tabout.file,
                             program         = "LTRharvest",
                             ltr.fasta       = NULL,
                             inner.seq.fasta = NULL,
                             data            = NULL,
                             similarity.bin  = 2,
                             min.sim         = NULL){
    
    if (!is.null(gff.file) & !is.null(data))
        stop ("Please only provide either a stored dataset or the path to the gff file that shall be imported.")
    
    if (!is.element(program, c("LTRharvest", "LTRdigest")))
        stop ("Please choose a prediction returned by either LTRharvest or LTRdigest.")
    
    X4 <- X5 <- X9 <- annotation <- ltr_similarity <- chromosome <- NULL
    
    if (program == "LTRharvest"){
        if (!is.null(ltr.fasta)){
            
            LTRretrotransposon.fasta <- Biostrings::readDNAStringSet(ltr.fasta,"fasta")
            LTRretrotransposon.NAMES <- LTRretrotransposon.fasta@ranges@NAMES
            RawIntervals <- do.call(rbind,sapply(LTRretrotransposon.NAMES, function(x) noquote(stringr::str_split(stringr::str_replace(stringr::str_replace(stringr::str_extract(x,"[?<=\\[].*?[?=\\]]"),"\\[",""),"\\]",""),","))))
            colnames(RawIntervals) <- c("start","end")
            ChrID <- sapply(rownames(RawIntervals), function(y) stringr::str_split(y, " \\(")[[1]][1])
            
            LTRretrotransposon.df <- data.frame(chromosome = ChrID, start = as.numeric(LTRretrotransposon.df[ , "start"]), end = as.numeric(LTRretrotransposon.df[ , "end"]))
            LTRretrotransposon.df <- dplyr::mutate(LTRretrotransposon.df, width = (end - start) + 1)
            
        }
        
        if (!is.null(inner.seq.fasta)){
            
            BetweenLTRSeq.fasta <- Biostrings::readDNAStringSet(inner.seq.fasta,"fasta")
            BetweenLTRSeq.NAMES <- BetweenLTRSeq.fasta@ranges@NAMES
            BetweenLTRRawIntervals <- do.call(rbind,sapply(BetweenLTRSeq.NAMES, function(x) noquote(stringr::str_split(stringr::str_replace(stringr::str_replace(stringr::str_extract(x,"[?<=\\[].*?[?=\\]]"),"\\[",""),"\\]",""),","))))
            colnames(BetweenLTRRawIntervals) <- c("start","end")
            ChrID <- sapply(rownames(BetweenLTRRawIntervals), function(y) stringr::str_split(y, " \\(")[[1]][1])
            
            BetweenLTRSeq.df <- data.frame(chromosome = ChrID, start = as.numeric(LTRretrotransposon.df[ , "start"]), end = as.numeric(LTRretrotransposon.df[ , "end"]))
            BetweenLTRSeq.df <- dplyr::mutate(LTRretrotransposon.df, width = (end - start) + 1)
            
        }
        
        if (!is.null(gff.file)){
            # check whether or not comment lines are included in the 
            # gtf file : '#!'
            CommentCheck <- readLines(gff.file,n = 1000, warn = FALSE)
            FoundCommentLines <- sapply(CommentCheck, function(x) stringr::str_detect(x, "#."))
            if (any(FoundCommentLines)){
                RemoveFirstNLines <- max(which(FoundCommentLines))
                
                # Extract header information
                #CommentCheck[2:RemoveFirstNLines]
                
            } else {
                RemoveFirstNLines <- 0
            }
            
            # read gff file content: comment = "###"
            AnnotationFile <- readr::read_tsv(gff.file, col_names = FALSE, skip = 0, comment = "#")
            cat("\n")
            cat("Input: ",gff.file," -> Row Number: ",nrow(AnnotationFile))
            cat("\n")
            AnnotationFile <- na.omit(AnnotationFile)
            cat("Remove 'NA' -> New Row Number: ",nrow(AnnotationFile))
            cat("\n")
        }
        
        
        if (!is.null(data))
            AnnotationFile <- data
        
        # determine the width of the predicted intervals
        AnnotationFile <- dplyr::mutate(AnnotationFile, width = (X5 - X4) + 1)
        colnames(AnnotationFile) <- c("chromosome","pred_tool","annotation","start","end","X6","strand","X8","X9","width")
        
        ### Post-Processing of repeat_region
        FilteredAnnotationFile.repeat_region <- dplyr::filter(AnnotationFile, annotation == "repeat_region")
        
        # Extract ID feature in column X9 (predicted by LTRharvest for repeat_region)
        RepeatRegionID <- vector("character",nrow(FilteredAnnotationFile.repeat_region))
        RepeatRegionID <- sapply(seq_len(nrow(FilteredAnnotationFile.repeat_region)), function(location) {
            
            unlist(stringr::str_split(FilteredAnnotationFile.repeat_region[location, "X9"],"="))[2]
            
        })
        
        FilteredAnnotationFile.repeat_region <- dplyr::mutate(FilteredAnnotationFile.repeat_region, repeat_region = RepeatRegionID)
        FilteredAnnotationFile.repeat_region <- dplyr::select(FilteredAnnotationFile.repeat_region,-X9)
        
        cat("(1/5) Filtering for repeat regions has been finished.")
        cat("\n")
        
        ### Post-Processing of LTR_retrotransposons
        FilteredAnnotationFile.LTR_retrotransposon <- dplyr::filter(AnnotationFile, annotation == "LTR_retrotransposon")
        # Extract ID, Parent, seq_number, and ltr_similarity features from column X9 (predicted by LTRharvest for LTR_retrotransposons)
        LTR_retrotransposonPredictionFeatures <- as.data.frame(t(sapply(seq_len(nrow(FilteredAnnotationFile.LTR_retrotransposon)), 
                                                                        function(location) unlist(lapply(stringr::str_split(unlist(stringr::str_split(FilteredAnnotationFile.LTR_retrotransposon[location ,"X9"],";")),"="), function(x) x[2])))), 
                                                               colClasses = c(rep("character",2),rep("numeric",2)))
        
        LTR_retrotransposonPredictionFeatures[ , 3] <- as.numeric(as.vector(LTR_retrotransposonPredictionFeatures[ , 3]))
        LTR_retrotransposonPredictionFeatures[ , 4] <- as.numeric(as.vector(LTR_retrotransposonPredictionFeatures[ , 4]))
        
        colnames(LTR_retrotransposonPredictionFeatures) <- c("ID", "repeat_region","ltr_similarity", "seq_number")
        LTR_retrotransposonPredictionFeatures <- cbind(FilteredAnnotationFile.LTR_retrotransposon[ , -9],LTR_retrotransposonPredictionFeatures)
        
        
        #     for (i in seq_len(length(table(LTR_retrotransposonPredictionFeatures[ , "chromosome"])))){
        #         if (i %in% 1:10)
        #             LTR_retrotransposonPredictionFeatures[which(LTR_retrotransposonPredictionFeatures[ , "chromosome"] == as.character(i)), "chromosome"] <- paste0("Chr",i) 
        #     }
        
        min.similarity <- vector("numeric",1)
        max.similarity <- vector("numeric",1)
        
        if (is.null(min.sim))
            min.similarity <- floor(min(LTR_retrotransposonPredictionFeatures[ , "ltr_similarity"]))
        if (!is.null(min.sim))
            min.similarity <- min.sim
        
        #     if ((min.similarity > min.sim) | (max.sim > max.similarity))
        #         stop ("Please enter a similarity interval between [0,100] that is present in the dataset!")
        #    evenCheck <- (100 - min.similarity) %% similarity.bin
        LTR_retrotransposonPredictionFeatures <- dplyr::mutate(LTR_retrotransposonPredictionFeatures, 
                                                               similarity = cut(ltr_similarity,
                                                                                rev(seq(100,min.similarity,-similarity.bin)),
                                                                                include.lowest = TRUE,
                                                                                right          = TRUE))
        
        #LTR_retrotransposonPredictionFeatures[ , "chromosome"] <- factor(LTR_retrotransposonPredictionFeatures[ , "chromosome"], ordered = TRUE, levels = paste0("Chr",1:10))
        
        FilteredAnnotationFile.repeat_region <- suppressWarnings(dplyr::right_join(FilteredAnnotationFile.repeat_region, LTR_retrotransposonPredictionFeatures[ , c("repeat_region","ltr_similarity")], by = "repeat_region"))
        FilteredAnnotationFile.repeat_region <- dplyr::mutate(FilteredAnnotationFile.repeat_region, 
                                                              similarity = cut(ltr_similarity,
                                                                               rev(seq(100,min.similarity,-similarity.bin)),
                                                                               include.lowest = TRUE,
                                                                               right          = TRUE))
        
        
        cat("(2/5) Filtering for LTR retrotransposons has been finished.")
        cat("\n")
        
        ### Post-Processing of inverted_repeat
        FilteredAnnotationFile.inverted_repeat <- dplyr::filter(AnnotationFile, annotation == "inverted_repeat")
        
        # Extract Parent feature in column X9 (predicted by LTRharvest for inverted_repeat)
        InvertedRepeatParentID <- vector("character",nrow(FilteredAnnotationFile.inverted_repeat))
        InvertedRepeatParentID <- sapply(seq_len(nrow(FilteredAnnotationFile.inverted_repeat)), function(location) {
            
            unlist(stringr::str_split(FilteredAnnotationFile.inverted_repeat[location, "X9"],"="))[2]
            
        })
        
        FilteredAnnotationFile.inverted_repeat <- dplyr::mutate(FilteredAnnotationFile.inverted_repeat, repeat_region = InvertedRepeatParentID)
        FilteredAnnotationFile.inverted_repeat <- dplyr::select(FilteredAnnotationFile.inverted_repeat,-X9)
        
        cat("(3/5) Filtering for inverted repeats has been finished.")
        cat("\n")
        
        ### Post-Processing of long_terminal_repeat
        FilteredAnnotationFile.long_terminal_repeat <- dplyr::filter(AnnotationFile, annotation == "long_terminal_repeat")
        # Extract Parent feature in column X9 (predicted by LTRharvest for long_terminal_repeat)
        LTRParentID <- vector("character",nrow(FilteredAnnotationFile.long_terminal_repeat))
        LTRParentID <- sapply(seq_len(nrow(FilteredAnnotationFile.long_terminal_repeat)), function(location) {
            
            unlist(stringr::str_split(FilteredAnnotationFile.long_terminal_repeat[location, "X9"],"="))[2]
            
        })
        
        FilteredAnnotationFile.long_terminal_repeat <- dplyr::mutate(FilteredAnnotationFile.long_terminal_repeat, ID = LTRParentID)
        FilteredAnnotationFile.long_terminal_repeat <- dplyr::select(FilteredAnnotationFile.long_terminal_repeat,-X9)
        
        FilteredAnnotationFile.long_terminal_repeat <- suppressWarnings(dplyr::right_join(FilteredAnnotationFile.long_terminal_repeat, LTR_retrotransposonPredictionFeatures[ , c("ID","ltr_similarity")], by = "ID"))
        FilteredAnnotationFile.long_terminal_repeat <- dplyr::mutate(FilteredAnnotationFile.long_terminal_repeat, 
                                                                     similarity = cut(ltr_similarity,
                                                                                      rev(seq(100,min.similarity,-similarity.bin)),
                                                                                      include.lowest = TRUE,
                                                                                      right          = TRUE))
        
        FilteredAnnotationFile.long_terminal_repeat$ID <- factor(FilteredAnnotationFile.long_terminal_repeat$ID, levels=unique(FilteredAnnotationFile.long_terminal_repeat$ID))
        FilteredAnnotationFile.long_terminal_repeat <- dplyr::mutate(FilteredAnnotationFile.long_terminal_repeat, ltr_order = do.call(rbind,lapply(split(FilteredAnnotationFile.long_terminal_repeat,FilteredAnnotationFile.long_terminal_repeat$ID), function(x) if (x[1, "end"] < x[2, "start"]) return (rbind("left","right")) else return(rbind("right","left")))))
        
        
        cat("(4/5) Filtering for LTRs has been finished.")
        cat("\n")
        
        ### Post-Processing of target_site_duplication
        FilteredAnnotationFile.target_site_duplication <- dplyr::filter(AnnotationFile, annotation == "target_site_duplication")
        # Extract Parent feature in column X9 (predicted by LTRharvest for target_site_duplication)
        TSDParentID <- vector("character",nrow(FilteredAnnotationFile.target_site_duplication))
        TSDParentID <- sapply(seq_len(nrow(FilteredAnnotationFile.target_site_duplication)), function(location) {
            
            unlist(stringr::str_split(FilteredAnnotationFile.target_site_duplication[location, "X9"],"="))[2]
            
        })
        
        FilteredAnnotationFile.target_site_duplication <- dplyr::mutate(FilteredAnnotationFile.target_site_duplication, repeat_region = TSDParentID)
        FilteredAnnotationFile.target_site_duplication <- dplyr::select(FilteredAnnotationFile.target_site_duplication,-X9)
        
        cat("(5/5) Filtering for target site duplication has been finished.")
        cat("\n")
        # rename the chromosome number
        FilteredAnnotationFile.repeat_region <- dplyr::mutate(FilteredAnnotationFile.repeat_region, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        LTR_retrotransposonPredictionFeatures <- dplyr::mutate(LTR_retrotransposonPredictionFeatures, seq_number =  chromosome)
        LTR_retrotransposonPredictionFeatures <- dplyr::mutate(LTR_retrotransposonPredictionFeatures, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        FilteredAnnotationFile.long_terminal_repeat <- dplyr::mutate(FilteredAnnotationFile.long_terminal_repeat, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        FilteredAnnotationFile.inverted_repeat <- dplyr::mutate(FilteredAnnotationFile.inverted_repeat, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        FilteredAnnotationFile.target_site_duplication <- dplyr::mutate(FilteredAnnotationFile.target_site_duplication, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        
        return(list( repeat.region           = FilteredAnnotationFile.repeat_region, 
                     ltr.retrotransposon     = LTR_retrotransposonPredictionFeatures,
                     ltr                     = FilteredAnnotationFile.long_terminal_repeat,
                     inverted_repeat         = FilteredAnnotationFile.inverted_repeat,
                     target.site.duplication = FilteredAnnotationFile.target_site_duplication) )
        
    }
    
    
    if (program == "LTRdigest"){
        
       # read gff file content: comment = "###"
       AnnotationFile <- readr::read_tsv(gff.file, col_names = FALSE, skip = 0, comment = "#")
       cat("\n")
       cat("Input: ",gff.file," -> Row Number: ",nrow(AnnotationFile))
       cat("\n")
       AnnotationFile <- na.omit(AnnotationFile)
       cat("Remove 'NA' -> New Row Number: ",nrow(AnnotationFile))
       cat("\n")
        
        
        if (!is.null(data))
            AnnotationFile <- data
        
        # determine the width of the predicted intervals
        AnnotationFile <- dplyr::mutate(AnnotationFile, width = (X5 - X4) + 1)
        colnames(AnnotationFile) <- c("chromosome","pred_tool","annotation","start","end","X6","strand","X8","X9","width")
        
        ### Post-Processing of repeat_region
        FilteredAnnotationFile.repeat_region <- dplyr::filter(AnnotationFile, annotation == "repeat_region")
        
        # Extract ID feature in column X9 (predicted by LTRharvest for repeat_region)
        RepeatRegionID <- vector("character",nrow(FilteredAnnotationFile.repeat_region))
        RepeatRegionID <- sapply(seq_len(nrow(FilteredAnnotationFile.repeat_region)), function(location) {
            
            unlist(stringr::str_split(FilteredAnnotationFile.repeat_region[location, "X9"],"="))[2]
            
        })
        
        FilteredAnnotationFile.repeat_region <- dplyr::mutate(FilteredAnnotationFile.repeat_region, repeat_region = RepeatRegionID)
        FilteredAnnotationFile.repeat_region <- dplyr::select(FilteredAnnotationFile.repeat_region,-X9)
        
        cat("(1/8) Filtering for repeat regions has been finished.")
        cat("\n")
        
        ### Post-Processing of LTR_retrotransposons
        FilteredAnnotationFile.LTR_retrotransposon <- dplyr::filter(AnnotationFile, annotation == "LTR_retrotransposon")
        # Extract ID, Parent, seq_number, and ltr_similarity features from column X9 (predicted by LTRharvest for LTR_retrotransposons)
        LTR_retrotransposonPredictionFeatures <- as.data.frame(t(sapply(seq_len(nrow(FilteredAnnotationFile.LTR_retrotransposon)), 
                                                                        function(location) unlist(lapply(stringr::str_split(unlist(stringr::str_split(FilteredAnnotationFile.LTR_retrotransposon[location ,"X9"],";")),"="), function(x) x[2])))), 
                                                               colClasses = c(rep("character",2),rep("numeric",2)))
        
        LTR_retrotransposonPredictionFeatures[ , 3] <- as.numeric(as.vector(LTR_retrotransposonPredictionFeatures[ , 3]))
        LTR_retrotransposonPredictionFeatures[ , 4] <- as.numeric(as.vector(LTR_retrotransposonPredictionFeatures[ , 4]))
        
        colnames(LTR_retrotransposonPredictionFeatures) <- c("ID", "repeat_region","ltr_similarity", "seq_number")
        LTR_retrotransposonPredictionFeatures <- cbind(FilteredAnnotationFile.LTR_retrotransposon[ , -9],LTR_retrotransposonPredictionFeatures)
        
        
        #     for (i in seq_len(length(table(LTR_retrotransposonPredictionFeatures[ , "chromosome"])))){
        #         if (i %in% 1:10)
        #             LTR_retrotransposonPredictionFeatures[which(LTR_retrotransposonPredictionFeatures[ , "chromosome"] == as.character(i)), "chromosome"] <- paste0("Chr",i) 
        #     }
        
        min.similarity <- vector("numeric",1)
        max.similarity <- vector("numeric",1)
        
        if (is.null(min.sim))
            min.similarity <- floor(min(LTR_retrotransposonPredictionFeatures[ , "ltr_similarity"]))
        if (!is.null(min.sim))
            min.similarity <- min.sim
        
        #     if ((min.similarity > min.sim) | (max.sim > max.similarity))
        #         stop ("Please enter a similarity interval between [0,100] that is present in the dataset!")
        #    evenCheck <- (100 - min.similarity) %% similarity.bin
        LTR_retrotransposonPredictionFeatures <- dplyr::mutate(LTR_retrotransposonPredictionFeatures, 
                                                               similarity = cut(ltr_similarity,
                                                                                rev(seq(100,min.similarity,-similarity.bin)),
                                                                                include.lowest = TRUE,
                                                                                right          = TRUE))
        
        #LTR_retrotransposonPredictionFeatures[ , "chromosome"] <- factor(LTR_retrotransposonPredictionFeatures[ , "chromosome"], ordered = TRUE, levels = paste0("Chr",1:10))
        
        FilteredAnnotationFile.repeat_region <- suppressWarnings(dplyr::right_join(FilteredAnnotationFile.repeat_region, LTR_retrotransposonPredictionFeatures[ , c("repeat_region","ltr_similarity")], by = "repeat_region"))
        FilteredAnnotationFile.repeat_region <- dplyr::mutate(FilteredAnnotationFile.repeat_region, 
                                                              similarity = cut(ltr_similarity,
                                                                               rev(seq(100,min.similarity,-similarity.bin)),
                                                                               include.lowest = TRUE,
                                                                               right          = TRUE))
        
        LTRdigest.tabout.file <- read.tabout(tabout.file)
        LTR_retrotransposonPredictionFeatures <- dplyr::inner_join(LTR_retrotransposonPredictionFeatures, LTRdigest.tabout.file, by = c("start" = "element_start", "end" = "element_end"))

        cat("(2/8) Filtering for LTR retrotransposons has been finished.")
        cat("\n")
        
        ### Post-Processing of inverted_repeat
        FilteredAnnotationFile.inverted_repeat <- dplyr::filter(AnnotationFile, annotation == "inverted_repeat")
        
        # Extract Parent feature in column X9 (predicted by LTRharvest for inverted_repeat)
        InvertedRepeatParentID <- vector("character",nrow(FilteredAnnotationFile.inverted_repeat))
        InvertedRepeatParentID <- sapply(seq_len(nrow(FilteredAnnotationFile.inverted_repeat)), function(location) {
            
            unlist(stringr::str_split(FilteredAnnotationFile.inverted_repeat[location, "X9"],"="))[2]
            
        })
        
        FilteredAnnotationFile.inverted_repeat <- dplyr::mutate(FilteredAnnotationFile.inverted_repeat, repeat_region = InvertedRepeatParentID)
        FilteredAnnotationFile.inverted_repeat <- dplyr::select(FilteredAnnotationFile.inverted_repeat,-X9)
        
        cat("(3/8) Filtering for inverted repeats has been finished.")
        cat("\n")
        
        ### Post-Processing of long_terminal_repeat
        FilteredAnnotationFile.long_terminal_repeat <- dplyr::filter(AnnotationFile, annotation == "long_terminal_repeat")
        # Extract Parent feature in column X9 (predicted by LTRharvest for long_terminal_repeat)
        LTRParentID <- vector("character",nrow(FilteredAnnotationFile.long_terminal_repeat))
        LTRParentID <- sapply(seq_len(nrow(FilteredAnnotationFile.long_terminal_repeat)), function(location) {
            
            unlist(stringr::str_split(FilteredAnnotationFile.long_terminal_repeat[location, "X9"],"="))[2]
            
        })
        
        FilteredAnnotationFile.long_terminal_repeat <- dplyr::mutate(FilteredAnnotationFile.long_terminal_repeat, ID = LTRParentID)
        FilteredAnnotationFile.long_terminal_repeat <- dplyr::select(FilteredAnnotationFile.long_terminal_repeat,-X9)
        
        FilteredAnnotationFile.long_terminal_repeat <- suppressWarnings(dplyr::right_join(FilteredAnnotationFile.long_terminal_repeat, LTR_retrotransposonPredictionFeatures[ , c("ID","ltr_similarity")], by = "ID"))
        FilteredAnnotationFile.long_terminal_repeat <- dplyr::mutate(FilteredAnnotationFile.long_terminal_repeat, 
                                                                     similarity = cut(ltr_similarity,
                                                                                      rev(seq(100,min.similarity,-similarity.bin)),
                                                                                      include.lowest = TRUE,
                                                                                      right          = TRUE))
        
        FilteredAnnotationFile.long_terminal_repeat$ID <- factor(FilteredAnnotationFile.long_terminal_repeat$ID, levels=unique(FilteredAnnotationFile.long_terminal_repeat$ID))
        FilteredAnnotationFile.long_terminal_repeat <- dplyr::mutate(FilteredAnnotationFile.long_terminal_repeat, ltr_order = do.call(rbind,lapply(split(FilteredAnnotationFile.long_terminal_repeat,FilteredAnnotationFile.long_terminal_repeat$ID), function(x) if (x[1, "end"] < x[2, "start"]) return (rbind("left","right")) else return(rbind("right","left")))))
        
        
        cat("(4/8) Filtering for LTRs has been finished.")
        cat("\n")
        
        ### Post-Processing of target_site_duplication
        FilteredAnnotationFile.target_site_duplication <- dplyr::filter(AnnotationFile, annotation == "target_site_duplication")
        # Extract Parent feature in column X9 (predicted by LTRharvest for target_site_duplication)
        TSDParentID <- vector("character",nrow(FilteredAnnotationFile.target_site_duplication))
        TSDParentID <- sapply(seq_len(nrow(FilteredAnnotationFile.target_site_duplication)), function(location) {
            
            unlist(stringr::str_split(FilteredAnnotationFile.target_site_duplication[location, "X9"],"="))[2]
            
        })
        
        FilteredAnnotationFile.target_site_duplication <- dplyr::mutate(FilteredAnnotationFile.target_site_duplication, repeat_region = TSDParentID)
        FilteredAnnotationFile.target_site_duplication <- dplyr::select(FilteredAnnotationFile.target_site_duplication,-X9)
        
        cat("(5/8) Filtering for target site duplication has been finished.")
        cat("\n")
        
        
        
        ### Post-Processing of primer_binding_site
        FilteredAnnotationFile.primer_binding_site <- dplyr::filter(AnnotationFile, annotation == "primer_binding_site")
        
        # Extract ID, Parent, seq_number, and ltr_similarity features from column X9 (predicted by LTRharvest for LTR_retrotransposons)
        PBSPredictionFeatures <- as.data.frame(t(sapply(seq_len(nrow(FilteredAnnotationFile.primer_binding_site)), 
                                                                        function(location) unlist(lapply(stringr::str_split(unlist(stringr::str_split(FilteredAnnotationFile.primer_binding_site[location ,"X9"],";")),"="), function(x) x[2])))), 
                                                               colClasses = c(rep("character",2),rep("numeric",3)))
        
        PBSPredictionFeatures[ , 3] <- as.numeric(as.vector(PBSPredictionFeatures[ , 3]))
        PBSPredictionFeatures[ , 4] <- as.numeric(as.vector(PBSPredictionFeatures[ , 4]))
        PBSPredictionFeatures[ , 5] <- as.numeric(as.vector(PBSPredictionFeatures[ , 5]))
        
        colnames(PBSPredictionFeatures) <- c("ID", "trna","trnaoffset", "pbsoffset","edist")
        FilteredAnnotationFile.primer_binding_site <- cbind(FilteredAnnotationFile.primer_binding_site[ , -9],PBSPredictionFeatures)
        
        FilteredAnnotationFile.primer_binding_site <- suppressWarnings(dplyr::right_join(FilteredAnnotationFile.primer_binding_site, LTR_retrotransposonPredictionFeatures[ , c("ID","ltr_similarity")], by = "ID"))
        FilteredAnnotationFile.primer_binding_site <- dplyr::mutate(FilteredAnnotationFile.primer_binding_site, 
                                                              similarity = cut(ltr_similarity,
                                                                               rev(seq(100,min.similarity,-similarity.bin)),
                                                                               include.lowest = TRUE,
                                                                               right          = TRUE))
        
        
        cat("(6/8) Filtering for primer binding site has been finished.")
        cat("\n")
        
        
        ### Post-Processing of protein_match
        FilteredAnnotationFile.protein_match <- dplyr::filter(AnnotationFile, annotation == "protein_match")
        
        # Extract ID, Parent, seq_number, and ltr_similarity features from column X9 (predicted by LTRharvest for LTR_retrotransposons)
        ProteinMatchPredictionFeatures <- as.data.frame(t(sapply(seq_len(nrow(FilteredAnnotationFile.protein_match)), 
                                                        function(location) unlist(lapply(stringr::str_split(unlist(stringr::str_split(FilteredAnnotationFile.protein_match[location ,"X9"],";")),"="), function(x) x[2])))), 
                                               colClasses = c("character","numeric","character"))
        
        ProteinMatchPredictionFeatures[ , 2] <- as.numeric(as.vector(ProteinMatchPredictionFeatures[ , 2]))
        
        
        colnames(ProteinMatchPredictionFeatures) <- c("ID", "reading_frame","name")
        FilteredAnnotationFile.protein_match <- cbind(FilteredAnnotationFile.protein_match[ , -9],ProteinMatchPredictionFeatures)
        
        FilteredAnnotationFile.protein_match <- suppressWarnings(dplyr::right_join(FilteredAnnotationFile.protein_match, LTR_retrotransposonPredictionFeatures[ , c("ID","ltr_similarity","width")], by = "ID"))
        FilteredAnnotationFile.protein_match <- dplyr::mutate(FilteredAnnotationFile.protein_match, 
                                                                    similarity = cut(ltr_similarity,
                                                                                     rev(seq(100,min.similarity,-similarity.bin)),
                                                                                     include.lowest = TRUE,
                                                                                     right          = TRUE))
        colnames(FilteredAnnotationFile.protein_match)[9] <- "match_width"
        colnames(FilteredAnnotationFile.protein_match)[14] <- "width"
        
        cat("(7/8) Filtering for primer binding site has been finished.")
        cat("\n")
        
        ### Post-Processing of RR_tract
        FilteredAnnotationFile.RR_tract <- dplyr::filter(AnnotationFile, annotation == "RR_tract")
        
        # Extract ID, Parent, seq_number, and ltr_similarity features from column X9 (predicted by LTRharvest for LTR_retrotransposons)
        RRtractPredictionFeatures <- as.data.frame(sapply(seq_len(nrow(FilteredAnnotationFile.RR_tract)), 
                                                                 function(location) unlist(stringr::str_split(FilteredAnnotationFile.RR_tract[location ,"X9"],"="))[2]), 
                                                        colClasses = "character")
        
        colnames(RRtractPredictionFeatures) <- "ID"
        FilteredAnnotationFile.RR_tract <- cbind(FilteredAnnotationFile.RR_tract[ , -9],RRtractPredictionFeatures)
        
        FilteredAnnotationFile.RR_tract <- suppressWarnings(dplyr::right_join(FilteredAnnotationFile.RR_tract, LTR_retrotransposonPredictionFeatures[ , c("ID","ltr_similarity")], by = "ID"))
        FilteredAnnotationFile.RR_tract <- dplyr::mutate(FilteredAnnotationFile.RR_tract, 
                                                              similarity = cut(ltr_similarity,
                                                                               rev(seq(100,min.similarity,-similarity.bin)),
                                                                               include.lowest = TRUE,
                                                                               right          = TRUE))
        
        cat("(8/8) Filtering for primer RR tract has been finished.")
        cat("\n")
        
        # rename the chromosome number
        FilteredAnnotationFile.repeat_region <- dplyr::mutate(FilteredAnnotationFile.repeat_region, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        LTR_retrotransposonPredictionFeatures <- dplyr::mutate(LTR_retrotransposonPredictionFeatures, seq_number =  chromosome)
        LTR_retrotransposonPredictionFeatures <- dplyr::mutate(LTR_retrotransposonPredictionFeatures, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        FilteredAnnotationFile.long_terminal_repeat <- dplyr::mutate(FilteredAnnotationFile.long_terminal_repeat, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        FilteredAnnotationFile.inverted_repeat <- dplyr::mutate(FilteredAnnotationFile.inverted_repeat, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        FilteredAnnotationFile.target_site_duplication <- dplyr::mutate(FilteredAnnotationFile.target_site_duplication, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        FilteredAnnotationFile.primer_binding_site <- dplyr::mutate(FilteredAnnotationFile.primer_binding_site, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        FilteredAnnotationFile.protein_match <- dplyr::mutate(FilteredAnnotationFile.protein_match, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        FilteredAnnotationFile.RR_tract <- dplyr::mutate(FilteredAnnotationFile.RR_tract, chromosome =  as.numeric(stringr::str_replace(chromosome,"seq","")) + 1)
        
        return(list( repeat.region           = FilteredAnnotationFile.repeat_region, 
                     ltr.retrotransposon     = LTR_retrotransposonPredictionFeatures,
                     ltr                     = FilteredAnnotationFile.long_terminal_repeat,
                     inverted_repeat         = FilteredAnnotationFile.inverted_repeat,
                     target.site.duplication = FilteredAnnotationFile.target_site_duplication,
                     pbs                     = FilteredAnnotationFile.primer_binding_site,
                     protein.match           = FilteredAnnotationFile.protein_match,
                     RR_tract                = FilteredAnnotationFile.RR_tract))
        
    }
}


