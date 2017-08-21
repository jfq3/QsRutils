#' Check Variance
#'
#' Tests for Heterogeneity of Variances in make_comparisons Result
#'
#' @param otu.pc.transformed An OTU matrix of transformed data.
#' @param group.vector A vector of treatments.
#'
#' @return Prints test results to the console.
#'
#' @export
#'
#' @seealso make_comparisons
#'
#' @examples
#'
check_var <- function(otu.pc.transformed, group.vector) {
  for (i in 1:nrow(otu.pc.transformed)) {
    rslt <- fligner.test(otu.pc.transformed[ i, ], g = group.vector)
    print((rownames(otu.pc.transformed))[i])
    print(rslt)
  }
}
