AlignmentMatrix <- function(LTRpred.seqs){
  # LTRpred predicted LTR transposon sequences
  #pred.seqs <- read.seqs(LTRpred.seqs)
  GlobalAlignmentScoreMatrix <- matrix(NA_real_, ncol = length(LTRpred.seqs), nrow = length(LTRpred.seqs))
  
  for (i in 1:length(LTRpred.seqs)){
    for(j in 1:length(LTRpred.seqs)){
      if (is.na(GlobalAlignmentScoreMatrix[i,j])){
        GlobalAlignmentScoreMatrix[i,j] <- NW.Align(LTRpred.seqs[i],LTRpred.seqs[j])@score
        GlobalAlignmentScoreMatrix[j,i] <- GlobalAlignmentScoreMatrix[i,j]
      }
    }
  }
  
#   colnames(GlobalAlignmentScoreMatrix) <- names(LTRpred.seqs)
#   rownames(GlobalAlignmentScoreMatrix) <- names(LTRpred.seqs)
  
  return (GlobalAlignmentScoreMatrix)
}