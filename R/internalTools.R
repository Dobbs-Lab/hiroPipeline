#Internal Tools
#These are helper functions for several user functions;
#they are not intended to be called directly by the user

#' readCifLines
#'
#' @param inFile Cif file to be read in
#'
#' @return cifContents The contents of the CIF file
#' @export
#'
#' @examples

readCifLines <- function(inFile){
  cifContents <- readLines(con      = inFile,
                           n        = -1L,
                           ok       = TRUE,
                           warn     = TRUE,
                           encoding = "unknown",
                           skipNul  = FALSE)
  return(cifContents)
}

#' stripWhiteSpace
#'
#' This function removes all white space from character vectors. If it is passed a data frame, it will remove all white space from all columns with character data types.
#'
#' @param wsco An object to be stripped of white space
#'
#' @return stripped The object, stripped of white space
#' @export
#'
#' @examples
#' stripWhiteSpace(c("red", " gre en ", "  ", " blue "))
#' V1 <- c("  1", " 2 ", "    3")
#' V2 <- c(1, 2, 3)
#' dummyDF <- data.frame(V1, V2, stringsAsFactors=FALSE)
#' dummyDF
#' stripWhiteSpace(dummyDF)

stripWhiteSpace <- function(wsco) UseMethod("stripWhiteSpace")

#Default method; will not run the method unless it is passed a data frame with at least one column of characters or a character vector
stripWhiteSpace.default <- function(wsco){
  stop("Object is not a character vector")
}

#For handling character vectors
stripWhiteSpace.character <- function(wsco){
  stripped <- gsub('\\s+', '', wsco)
  return(stripped)
}

#For handling data frames with at least one character column type
stripWhiteSpace.data.frame <- function(wsco){
  count <- 0
  dDF <- wsco
  for(i in 1:ncol(wsco)){
    if(class(wsco[,i])=="character"){
      dDF[,i] <- stripWhiteSpace(wsco[,i])
      count <- count + 1
    }
  }
  #Determine if any of the columns are character vectors; stop execution if they are not
  if(count == 0){
    stop("Data frame has no character columns")
  }
  return(dDF)
}



#' cleanNonFastaInput
#'
#' @param sequence A sequence (such as from ENSEMBL or ENCODE) copy/pasted
#' from a web page with indexing lines and spaces
#'
#' @return The ACGT characters in sequence, with whitespace and numeric characters removed
#' @export
#'
#' @examples
#' sequence <- "1  ACGTACGTACGTACGTACGT 20
#'              21 ACGTACGTACGTACGTACGT 30"
#' cleanNonFastaInput(sequence)

cleanNonFastaInput <- function(sequence) UseMethod("cleanNonFastaInput")

cleanNonFastaInput.default <- function(sequence){
  stop("Input is not a character vector.")
}


cleanNonFastaInput.character <- function(sequence){
  #Remove whitespace
  cleaned <- stripWhiteSpace(sequence)
  #Remove numbers
  cleaned <- gsub('[0-9]', '', cleaned)
  return(cleaned)
}




#' complement
#'
#' This function takes a DNA or RNA sequence as input (along with a parameter specifying the type of sequence) and outputs the complement of the input sequence. E.g., "ATTG" will return "TAAC" if type = "DNA" and "UAAC" if type = "RNA"
#'
#' @param seq A DNA or RNA sequence from which to generate a complement string
#' @param type Default is "DNA"; a DNA sequence can only contain "A", "C", "G", or "T" for the purposes of complement(). The other option is "RNA"; an RNA sequence can only contain "A", "C", "G", or "U" for the purposes of complement().
#'
#' @return compSeq The complement of the input sequence
#' @export
#'
#' @examples
#' seq <- "AAAATGGCGAAG"
#' type <- "DNA"
#' complement(seq, type)


#Define the complement function; the type of function used is dependent on the type of object seq is
complement <- function(seq, type){
  UseMethod("complement", seq)
}

#Default complement function kicks back an error if seq is not a character sequence or list of
#character sequences to prevent attempts to complement objects that are not sequences
complement.default <- function(seq, type){
  stop("Error: Cannot complement a non-character vector or list object. Please check input sequence.")
}

#Complement function for character sequences
complement.character <- function(seq, type){
  compSeq <- seq

  #Create a pattern to check the submitted sequence for non-nucleotide characters

  #heck to ensure that the submitted sequence has characters allowed in DNA and RNA sequences
  notAllowed <- "[^ACGTUacgtuNn-]"

  #Check if there are any non-nucleotide characters present in the sequence
  if(grepl(notAllowed, compSeq)){
    stop("Error: Sequence contains non-DNA/non-RNA characters. Please check input sequence; allowed characters are 'A', 'C', 'T', and 'G' for [type = 'DNA'], and 'A', 'C', 'G', and 'U' for [type = 'RNA'].")
  } else {
    if(type == "DNA"){
      if(grepl("U", compSeq)){
        stop("Error: Uracil detected in sequence. Please use [type = 'RNA'] or check input sequence.")
      } else {
        fromVal <- c("A", "C", "G", "T", "a", "c", "g", "t", "N", "n", "-")
        toVal   <- c("T", "G", "C", "A", "t", "g", "c", "a", "N", "n", "-")
      }
    } else if(type == "RNA"){
      if(grepl("T", compSeq)){
        stop("Error: Thymine detected in sequence. Please use [type = 'DNA'] or check input sequence.")
      } else {
        fromVal <- c("A", "C", "G", "U", "a", "c", "g", "u", "N", "n", "-")
        toVal   <- c("U", "G", "C", "A", "u", "g", "c", "a", "N", "n", "-")
      }
    } else {
      stop("Error: Sequence type must be 'DNA' or 'RNA'.")
    }
  }

  #Replace the nucleotides with their complements
  compSeq <- plyr::mapvalues(unlist(strsplit(compSeq, split = "")),
                             from = fromVal,
                             to   = toVal,
                             warn_missing = FALSE)
  compSeq <- paste(compSeq, collapse = "")
  return(compSeq)
}


#Defines function to handle lists
complement.list <- function(seq, type){

  retList <- lapply(seq, complement, type = type)
  return(retList)
}



#' reverse
#'
#' This function takes a string as input and reverses the order of the string
#'
#' @param seq A string to reverse
#'
#' @return revSeq The seq string in reverse
#'
#' @examples
#' reverse(123456)
#' reverse("hello world")
#'
#' @export


#Define reverse function; the type of "seq" determines which version of "reverse()" is used
reverse <- function(seq){
  UseMethod("reverse", seq)
}

#Define default reverse function to throw an error fo seq objects that are not character strings or integers
reverse.default <- function(seq){
  stop("Error: Cannot reverse objects that are not character strings or integers. Please check input sequence.")
}

#Define reverse function for character strings and character vectors
reverse.character <- function(seq){
  #Determine if seq is a character vector (e.g., c("Hello", "World")) or a character string (e.g., "Hello")
  if(length(seq) > 1){
    #Apply reverse to each item on the list
    revSeq <- unlist(lapply(seq, reverse))
  } else {
    #Create a character string the same length as seq to restore the reversed sequence
    revSeq <- seq

    for(i in 1:nchar(seq)){
      #Overwrite the last character in revSeq with the first in seq, the second to last in revSeq with the second in seq, etc.
      curL <- substr(seq, i, i)
      substr(revSeq, (nchar(seq) + 1 - i), (nchar(seq) + 1 - i)) <- curL
    }
  }
  return(revSeq)
}

#Define reverse function for numeric sequences
reverse.numeric <- function(seq){
  #Convert the integer to a character sequence
  charSeq <- as.character(seq)

  #Call reverse() on the character sequence, and turn it back into a number afterwards
  revSeq <- as.numeric(reverse.character(charSeq))
  return(revSeq)
}

#Define reverse function for a list of sequences that need to be reversed
reverse.list <- function(seqList){

  #Apply reverse to each item on the list
  revSeqList <- lapply(seqList, reverse)

  return(revSeqList)
}


#' reverseComplement
#'
#' This function takes a DNA or RNA sequence as input and outputs the reverse complement of the sequence.
#'
#' @param seq A character vector from which to generate a reverse complement.
#' @param type Default is "DNA"; allowed characters are "A", "C", "G", and "T" (case insensitive). Other option is "RNA"; allowed characters are "A", "C", "G", and "U" (case insensitive.)
#'
#' @return seqRevComp The reverse complement of the input sequence
#' @export
#'
#' @examples
#' dnaSeq <- "AATGCC"
#' reverseComplement(dnaSeq)
#' rnaSeq <- "UUAGCC"
#' reverseComplement(rnaSeq, type = "RNA")


#Define reverseComplement function; set default type to "DNA"
reverseComplement <- function(seq, type = "DNA"){
  UseMethod("reverseComplement", seq)
}

#Prevent running reverseComplement on things that are not DNA/RNA sequences
reverseComplement.default <- function(seq, type = "DNA"){
  stop("Error: Input sequence is not a character vector. Please check input sequence.")
}

#Defines handling of character sequences
reverseComplement.character <- function(seq, type = "DNA"){
  #Get the reverse complement sequence
  seqRevComp <- complement(reverse(seq), type)
  return(seqRevComp)
}

#Define handling of lists
reverseComplement.list <- function(seq, type = "DNA"){
  retList <- lapply(seq, reverseComplement, type = type)

  return(unlist(retList))
}


