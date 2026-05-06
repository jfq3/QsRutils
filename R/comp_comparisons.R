#' Make Comparisons
#'
#' Calculates the treatment comparison portion of a table comparing relative
#' abundances of each taxon among treatments.
#'
#' @param otu.pc An OTU table of percentages.
#' @param otu.pc.trans An OTU table of transfromed data.
#' @param grps A vector of treatemnt groups for which to make comparisons.
#' @param p.adjust.method Adjustment method for multiple comparisons.
#' @param pool.sd A logical, whether or not to pool standard deviations.
#'
#' @return A data frame of differences in relative abundances among treatments.
#' @detail The data frame returned has taxa as row names. The first three column names are
#'   mean (relative abundance of the taxa), sd and F_value for comparisons among groups. The remaining column names are the groups. The group columns show the mean relative abundance for the group plus/minus the standard deviation and a compact letter display for the group.
#' @export
#' @importFrom stats pairwise.t.test
#' @importFrom stats sd
#' @examples
#' {
#' data("its.root")
#' temp1 <- comp_prepare_phyloseq(its.root)
#' temp2 <- comp_prepare_otu_table(temp1$expt.taxon.pc,
#'                                grps = "Label",
#'                                transformation = "sqrt_arc_sine")
#' comp_comparisons(otu.pc = temp2$otu.pc,
#'                  otu.pc.trans = temp2$otu.pc.trans,
#'                  grps = temp2$groups,
#'                  p.adjust.method = "BH",
#'                  pool.sd = TRUE)
#' }
#' 
comp_comparisons <- function(otu.pc, otu.pc.trans, grps, p.adjust.method = "BH", pool.sd = FALSE) {
  # Prepare comparison results matrix
  comp.sum <- matrix(NA, nrow=nrow(otu.pc.trans), ncol=(length(levels(grps)) + 1))
  colnames(comp.sum) <- c("Taxon", levels(grps))

  # Fill comparison result matrix
  for (i in 1:nrow(otu.pc.trans)) {
    taxrank.name <- as.character(rownames(otu.pc)[i]) # gets taxon name

    # Calculate means and standard deviations by row
    df <- data.frame(otu.pc[i, ], grps)
    colnames(df) <- c("value", "grps")
    avg <- tapply(df[ , 1], df[ , 2], mean) # calculates means
    std.dev <- tapply(df[ , 1], df[ , 2], sd) # calculates sd
    avg <- round(avg, 2)
    std.dev <- round(std.dev, 2)
    l <- paste(avg, "+/-", std.dev, sep = "")

    # Pairwise t-tests
    x <- stats::pairwise.t.test(otu.pc.trans[ i, ], g = grps, p.adjust.method = p.adjust.method, pool.sd = pool.sd)

    # Get letter assignments for results; allow for NAs, NaNs in results
    x <- get_groups(x)$p.matrix
    x <- replace(x, is.na(x), 0)
    ltrs <- multcompLetters(x)
    ltrs <- unname(ltrs$Letters)

    # Paste results together
    comp <- paste(l, ltrs, sep="")
    comp.sum[ i, ] <- c(taxrank.name, comp)
  }

  # Make first column rownames, then drop first column
  rownames(comp.sum) <- comp.sum[ , 1]
  comp.sum <- comp.sum[ , -1] |> 
    as.data.frame()

  return(comp.sum)
}
