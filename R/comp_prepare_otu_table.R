#' Prepare OTU Table
#'
#' Make OTU tables for making comparisons of relative abundances among treatments.
#'
#' @param expt.taxon.pc Phyloseq object from comp_prepare_phyloseq with percentages in the otu_table.
#' @param grps Factor in sample data for which to make comparisons.
#' @param transformation Transformation function to use.
#'
#' @details transformation may be "none" or a user-supplied function name in quotation marks or any of the built-it transformations("arc_sine", "log_arc_sine", or "sqrt_arc_sine"). The +sqrt_arc_sine" has generaally proven most effective.
#'
#' @return A list consisting of an OTU table with percentages, an OTU table with transformed data, and a vector of treatment groups.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' See vignette
#' }
#'
comp_prepare_otu_table <- function(expt.taxon.pc, grps = "Treatment", transformation = "sqrt_arc_sine"){
  # Extract percentage otu table, make rownames the taxa names.
  otu.pc <- veganotu(expt.taxon.pc)
  # taxrank <- rank_names(expt.pc)
  colnames(otu.pc) <- tax_table(expt.taxon.pc)

  # Apply transformation
  if (!(transformation == "none")) {
    otu.pc.trans <- apply(otu.pc, 1, transformation)
  } else {
    otu.pc.trans <- t(otu.pc)
  }

  otu.pc <- t(otu.pc)

  # Extract sample data table, get vector of groups
  sam <- vegansam(expt.taxon.pc)
  grps <- sam[ , grps]
  grps <- factor(grps, ordered = FALSE)

  return(list(otu.pc=otu.pc, otu.pc.trans=otu.pc.trans, groups=grps))
}
