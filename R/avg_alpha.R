#' Average Alpha Diversity
#'
#' Calculates alpha-diversity metrics from n samplings of an OTU table to a constant number of counts per sample.
#'
#' @param otu An OTU table as a data frame with samples as rows and taxa in columns.
#' @param sampling_depth The number of counts per sample in the sampled OTU table
#' @param replications The number of times the OTU table should be sampled.
#' @param sum_method Method (median or mean) for summarizing replication results.
#' 
#' @details
#' This function generates n rarefied community data frames (OTU tables) of a 
#' given sample size, calculates alpha-diversity metrics for each, and then 
#' summarizes the results.
#' The OTU data frame supplied must be in typical vegan format, that is with samples
#' as row names and taxa as column names.
#' The minimum row sum must be greater than the sampling depth.
#' 
#' @importFrom vegan diversity
#' @importFrom vegan rrarefy
#' @importFrom vegan specnumber
#' 
#' @seealso [vegann::rrarefy()]
#' 
#' @return Returns a dataframe with Shannon, Observed, Pielou, Simpson, InvSimpson for each sample in an OTU table.
#' @export
#'
#' @examples
#' {
#' data(BCI)
#' dim(BCI)
#' summary(rowSums(BCI))
#' otu <- BCI[rowSums(BCI)>400, ]
#' dim(otu)
#' summary(rowSums(otu))
#' avg_alpha(otu, sampling_depth = 400)
#' }
#'

#'
avg_alpha <- function(otu, sampling_depth, replications = 100, sum_method = "median") {
  a <- nrow(otu)
  b <- 5
  c <- replications
  a_names <- rownames(otu)
  b_names <- c("Shannon", "Observed", "Pielou", "Simpson", "InvSimpson")
  c_names <- paste0("rep", seq_along(1:replications))
  d <- array(0,
             dim = c(a, b, c),
             dimnames = list(a_names, b_names, c_names))
  for (i in seq(1:replications)) {
    otu_r <- vegan::rrarefy(otu, sampling_depth)
    d[, 1, i] <- vegan::diversity(otu_r, index = "shannon")
    d[, 2, i] <- vegan::specnumber(otu_r)
    d[, 3, i] <- d[, 1, i] / log(d[,2, i])
    d[, 4, i] <- vegan::diversity(otu_r, index = "simpson")
    d[, 5, i] <- vegan::diversity(otu_r, index = "invsimpson")
  }
  rslt <- apply(d, c(1,2), sum_method)
  return(rslt)
}
