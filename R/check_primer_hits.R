#' Check Primer Hits
#' 
#' Determine hits of all orientations of primers to paired sequence files.
#' 
#' @param path Path to paired sequence files
#' @param fwd_pattern Portion of file name that distinguishes forward read files. The default is "_R1.fastq"
#' @param rev_pattern Portion of the file name that distinguishes the reverse file. The default is "_R2.fastq"
#' @param fwd_primer Nucleotide sequence of the forward primer
#' @param rev_primer Nucleotide sequence of the reverse primet
#' @return A table of hits to the sequences by all primer orientations
#' @export
#' @details This function is for checking the effectiveness of primer trimming of ITS sequences.
#' Because the ITS region varies in length, it is possible for forward sequences to extend past the reverse primenr region and vice-versa, leadsing to what Robert Edgar calls staggered pairs fi the sequences are merged.
#' @details The fwd_pattern and rev_pattern must contain the file extension. The defaults are "_R1.fastq" and "_R2.fastq".
#' @details Default primers are ITS5 (forward) and ITS2 (reverse) from White et.al 1990.
#' @details ITS5: "GGAAGTAAAAGTCGTAACAAGG"
#' @details ITS2: "GCTGCGTTCTTCATCGATGC"
#' @importFrom Biostrings DNAString
#' @importFrom Biostrings complement
#' @importFrom Biostrings reverse
#' @importFrom Biostrings reverseComplement
#' @importFrom Biostrings vcountPattern
#' @importFrom dada2 filterAndTrim
#' @importFrom ShortRead readFastq
#' @importFrom ShortRead sread
#' @examples
#' # This takes a while
#' path <- system.file("extdata", package = "QsRutils")
#' check_primer_hits(
#' path = path,
#' fwd_pattern = "raw_1.fastq.gz",
#'   rev_pattern = "raw_2.fastq.gz",
#' fwd_primer = "GGAAGTAAAAGTCGTAACAAGG",
#' rev_primer = "GCTGCGTTCTTCATCGATGC"
#' )

check_primer_hits <- function(path = getwd(),
                              fwd_pattern= "_R1.fastq",
                              rev_pattern= "_R2.fastq",
                              fwd_primer = "GGAAGTAAAAGTCGTAACAAGG", 
                              rev_primer = "GCTGCGTTCTTCATCGATGC")
{
  
  fnFs <- sort(list.files(path, pattern = fwd_pattern, full.names = TRUE))
  fnRs <- sort(list.files(path, pattern = rev_pattern, full.names = TRUE))
  allOrients <- function(primer) {
    # Create all orientations of the input sequence
    dna <- Biostrings::DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, 
                 Complement = Biostrings::complement(dna), 
                 Reverse = Biostrings::reverse(dna),
                 RevComp = Biostrings::reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
  }
  
  FWD.orients <- allOrients(fwd_primer)
  REV.orients <- allOrients(rev_primer)
  
  fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
  fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
  
  dada2::filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)
  
  primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- Biostrings::vcountPattern(primer, sread(ShortRead::readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
  rslt <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]),
                FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]),
                REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]),
                REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
  return(rslt)
}
