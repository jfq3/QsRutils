#' Subset Distance Matrix
#'
#' Subsets a distance matrix.
#'
#' @param physeq An experiment level phyloseq object.
#' @param d.matrix A distance matrix.
#'
#' @return A distance matrix of smaller dimensions.
#' @export
#'
#' @details Some distance matrices take a long time to calculate for large data sets. This is especially true of unifrac and generalized unifrac distances calculated by GUniFracs. If distances are first calculated from data in a large experiment level phyloseq object and then it is desired to perform PERMANOVA (with adonis) on a subset of that object, this function provides a means of sub-setting the distance matrix so that it does not have to be calculated again for the subset data. The arguments are the distance matrix for the original phyloseq object and the smaller phyloseq object subset from the original.
#'
#' @references Chen J, Bittinger K, Charlson ES et al. (2012) Associating microbiome composition with environmental covariates using generalized UniFrac distances. Bioinformatics, 28, 2106-2113.
#'
#' @examples
#'
subset_dist <- function(physeq, d.matrix) {
  d.matrix <- as.matrix(d.matrix)
  keep <- sample_names(physeq)
  d.sub <- d.matrix[keep, keep]
  d.sub <- as.dist(d.sub)
  return(d.sub)
}
