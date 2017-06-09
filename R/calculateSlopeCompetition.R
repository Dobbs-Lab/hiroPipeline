#' calculateSlopeCompetition
#'
#' @param inData A table of two columns - first is geneID, second is wildtype sequence
#' @param cutSite The expected index of the nucleotide to the left of the cut; -1 if cut site is in middle of wildtype sequence
#' @param weight A weight factor for caculating pattern score; default is 20, as in Bae paper
#' @param top The number of entries to consider when generating the linear regression
#'
#' @return


#For each sequence in the list, calculate the slope competition
calculateSlopeCompetition <- function(inData, cutSite = -1, weight = 20, top = 10){
  #Create data frame to hold information about targets
  targetDF <- data.frame(geneId             = as.character(),
                         seq                = as.character(),
                         microhomologyScore = as.numeric(),
                         outOfFrameScore    = as.numeric(),
                         slopeMH2Plus       = as.numeric(),
                         slopeMH3Plus       = as.numeric(),
                         stringsAsFactors = FALSE)

  lengthWeight <- weight

  for(a in 1:nrow(inData)){
    if(cutSite == -1){
      #If the cut site isn't specified, assume it is in the middle of the input
      cutSite <- nchar(inData[a, 2])/2
    }

    #Get the pattern scores for every deletion pattern in the wildtype sequence
    patternScoreDF <- calculateBae(inData[a, 2], cutSite, lengthWeight)
    patternScoreDF <- patternScoreDF[order(-patternScoreDF$patternScore),]

    #Subset the out of frame scores
    outOfFrameInst <- patternScoreDF[which(patternScoreDF$delLength %% 3 != 0),]

    #Calculate the microhomology score according to Bae algorithm
    mhScore <- sum((nchar(patternScoreDF$microhomology) +
                      stringr:::str_count(toupper(patternScoreDF$microhomology), "G") +
                      stringr:::str_count(toupper(patternScoreDF$microhomology), "C")) *
                     (1/exp((patternScoreDF$delLength)/lengthWeight)) * 100)

    #Calculate the out-of-frame score according to Bae algorithm
    outOfFrameScore <- (sum((nchar(outOfFrameInst$microhomology) +
                               stringr:::str_count(toupper(outOfFrameInst$microhomology), "G") +
                               stringr:::str_count(toupper(outOfFrameInst$microhomology), "C")) *
                              (1/exp((outOfFrameInst$delLength)/lengthWeight)) * 100) / mhScore) * 100

    threePlus <- patternScoreDF[which(nchar(patternScoreDF$microhomology) >= 3),]
    if(nrow(threePlus) < top){
      end3 <- nrow(threePlus)
    } else {
      end3 <- top
    }
    top103Plus <- threePlus[1:end3,]

    if(nrow(patternScoreDF) < top){
      end2 <- nrow(patternScoreDF)
    } else {
      end2 <- top
    }

    #Perform linear regression on the top 10
    linModel3 <- lm(top103Plus$patternScore ~ seq(1:end3), data = top103Plus)
    linModel2 <- lm(patternScoreDF$patternScore[1:end2] ~ seq(1:end2), data = patternScoreDF)
    #Create data frame to hold information
    tempFrame <- data.frame(geneId = inData[a,1],
                            seq    = inData[a,2],
                            microhomologyScore = mhScore,
                            outOfFrameScore = outOfFrameScore,
                            slopeMH2Plus = linModel2$coefficients[2],
                            slopeMH3Plus = linModel3$coefficients[2],
                            stringsAsFactors = FALSE)

    targetDF <- rbind(targetDF, tempFrame)
  }

  return(targetDF)
}
