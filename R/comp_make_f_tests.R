#' Make F Tests
#'
#' Calclates omnibus F tests to be included in a table comparing relative abundances of each taxon among treatments.
#'
#' @param otu.pc.trans An OTU table of transformed data from comp_prepare_otu_table.
#' @param grps A vector of treatemnt groups for which to make comparisons.
#' @param var.equal Logical, whether or not to assume variances equal.
#'
#' @return A data frame of the F-test results.
#' @export
#'
#' @examples
#' \dontrun{
#' See vignette
#' }
#'
comp_make_f_tests <- function(otu.pc.trans, grps, var.equal = FALSE) {
  # Make global F tests

  # Prepare result data frame
  rslt <- as.data.frame(matrix(nrow=nrow(otu.pc.trans), ncol=3))
  colnames(rslt) <- c("Taxon", "F", "Prob")
  rslt[ , 1] <- rownames(otu.pc.trans)

  # Fill results data frame
  for (i in 1:nrow(otu.pc.trans)) {
    temp2 <- oneway.test(otu.pc.trans[ i, ] ~ grps, var.equal = var.equal)
    rslt[i,2] <- temp2$statistic
    rslt[i,3] <- temp2$p.value
  }

  rslt$sig <- lapply(rslt$Prob, asterix)
  rslt$F_value <- paste(round(rslt$F, 2), rslt$sig, sep="")
  rownames(rslt) <- rslt[ , 1]
  f.test.rslt <- rslt[ , -1]

  return(f.test.rslt)

}

