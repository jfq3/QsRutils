#' Check Variance
#'
#' Tests for Heterogeneity of Variances in make_comparisons Result
#'
#' @param otu.pc.transformed An OTU matrix of transformed data. Taxa are rows.
#' @param group.vector A vector of treatments.
#'
#' @return Prints test results to the console.
#'
#' @export
#'
#' @importFrom stats fligner.test
#' @seealso make_comparisons
#'
#' @examples
#' # Transform species matrix to proportion;
#' # Check variances for the first three species
#' # in the dune data set grouped by Management.
#' data(dune)
#' data(dune.env)
#' dune <- vegan::decostand(dune, method = "total")
#' dune <- dune[, 1:3]
#' dune <- t(dune)
#' check_var(dune, dune.env$Management)
#'
check_var <- function(otu.pc.transformed, group.vector) {
  for (i in 1:nrow(otu.pc.transformed)) {
    rslt <- stats::fligner.test(otu.pc.transformed[ i, ], g = group.vector)
    print((rownames(otu.pc.transformed))[i])
    print(rslt)
  }
}
