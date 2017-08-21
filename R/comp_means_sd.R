#' Calculate Means and Standard Deviations
#'
#' Calculates means and standard deviation for each taxon to be included in a table comparing relative abundances of each taxon among treatments.
#'
#' @param otu.pc An OTU table with data as percentages.
#'
#' @details The OTU table should be created with comp_prepare_otu_table.
#'
#' @return A data frame with means and standard deviations by taxon.
#' @export
#'
#' @examples
#'
comp_means_sd <- function(otu.pc) {
  # Calculate mean and standard deviation for each row of data.
  otu.pc.means <- apply(otu.pc, 1, mean)
  otu.pc.sd <- apply(otu.pc, 1, sd)
  otu.pc.sum <- data.frame(otu.pc.means, otu.pc.sd)
  colnames(otu.pc.sum) <- c("mean", "sd")
  otu.pc.sum <- round(otu.pc.sum, 2)
  return(otu.pc.sum)
}
