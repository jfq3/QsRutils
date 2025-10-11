#' @title Hash DNA sequences
#' @name hash_dna_seqs
#' @description Renames sequences in a multi-fasta file with the MD5 hash of each sequence.
#' @param seqs A list of DNA sequences
#' 
#' @returns DNA sequences renamed with their hashes.
#' @details Taxa names for the ASV table returned by the R version of DADA2 are the sequences themselves. This function renames them with the MD5 hash of each sequences so that the results are directly comparable QIIME2/DADA2 results.
#' @export
#' @examples
#' seqs <- "AGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAG"
#' hash_dna_seqs(seqs)
#' 
hash_dna_seqs <- function(seqs) {
  rslt <- insect::hash(seqs)
  return(rslt)
}