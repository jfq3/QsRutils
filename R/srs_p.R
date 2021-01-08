#' srs_p
#' 
#' Normalize sample counts using scaling with ranked subsampling (SRS)
#' 
#' @usage srs_p(p)
#' 
#' @param p a phyloseq object containing an OTU table
#'
#' @return a phyloseq object including an OTU table, all sample sums equal.
#' 
#' @export
#' 
#' @details This is an alternative to "rarefying" an OTU table to a constant sample size.  The phyloseq object submitted must be pruned to the desired sample size before uisng this function.
#' 
#' @references Beule L, Karlovsky P. Improved normalization of species count data in ecology by scaling with ranked subsampling (SRS): application to microbial communities. PeerJ. 2020;8:e9593.
#' 
srs_p <- function(p) {
  otu <- veganotu(p)
  taxa.names <- colnames(otu)
  otu <- t(otu)
  otu <- as.data.frame(otu)
  cmin <- min(sample_sums(p))
  otu.srs <- SRS(data=otu, Cmin = cmin)
  rownames(otu.srs) <- taxa.names
  new.otu <- otu_table(otu.srs, taxa_are_rows = TRUE)
  otu_table(p) <- new.otu
  return(p)
}