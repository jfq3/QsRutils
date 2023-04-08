#' Filter OTUs by Abundance
#'
#' Allows subsetting of a phyloseq object according to the relative abundance of OTUs in a minimal number of samples. Returns a logical vector of OTUs that are at least n\% of the sequences in at leas m samples.
#'
#' @param x A phyloseq object.
#' @param n Minimum percentage to keep OTU.
#' @param m Minimum number of samples.
#'
#' @return A logical vector of OTUs to keep.
#' @export
#'
#' @details The functions creates a logical vector to be used in subsetting a phyloseq object according to the relative abundance of OTUs in a given number of samples. For example, if n = 1 and m = 2, then the OTUs to be kept must represent at least 1\% of the sequences in at least 2 samples. The vector is then used as an argument to the phyloseq object `prune_taxa`.
#'
#' @examples
#' \dontrun{
#' prop_filter(expt, 1, 5)
#' }
#'
#' @importFrom phyloseq taxa_are_rows
#' @importFrom phyloseq otu_table
#' @importFrom vegan decostand
#'
prop_filter <- function(x, n, m) {
  test <- taxa_are_rows(x)
  if (test) {
    otu <- t(as(otu_table(x), "matrix"))
  }
  else {
    otu <- as(otu_table(x), "matrix")
  }
  otu <- 100*decostand(otu, "total")
  keep <- otu > n
  keep <- colSums(keep)>=m
  names(keep) <- colnames(otu)
  return(keep)
}
