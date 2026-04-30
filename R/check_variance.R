#' Check Variance
#'
#' Tests for Heterogeneity of Variances in make_comparisons Result
#'
#' @param otu.pc.transformed An OTU matrix of transformed data. Taxa are rows.
#' @param group.vector A vector of treatments.
#'
#' @return A data frame of class \code{"check_var"} with one row per taxon and
#'   columns \code{taxon}, \code{statistic}, \code{df}, and \code{p.value} from
#'   the Fligner-Killeen test. The object can be printed with \code{print()}.
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
#' data(dune, package = "vegan")
#' data(dune.env, package = "vegan")
#' dune <- vegan::decostand(dune, method = "total")
#' dune <- dune[, 1:3]
#' dune <- t(dune)
#' check_var(dune, dune.env$Management)
#'
check_var <- function(otu.pc.transformed, group.vector) {
  results <- lapply(seq_len(nrow(otu.pc.transformed)), function(i) {
    rslt <- stats::fligner.test(otu.pc.transformed[i, ], g = group.vector)
    data.frame(
      taxon     = rownames(otu.pc.transformed)[i],
      statistic = unname(rslt$statistic),
      df        = unname(rslt$parameter),
      p.value   = rslt$p.value,
      row.names = NULL
    )
  })
  out <- do.call(rbind, results)
  class(out) <- c("check_var", "data.frame")
  out
}

#' @export
#' @rdname check_var
print.check_var <- function(x, ...) {
  cat("Fligner-Killeen Test of Homogeneity of Variances\n\n")
  print(data.frame(x), ...)
  invisible(x)
}
