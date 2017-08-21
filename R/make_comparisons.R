#' Make Multiple Comparisons on Transformed Data
#'
#' Makes multiple comparisons of the relative abundances of taxa between treatment groups using the pairwise.t.test. Data may be transformed by a user supplied function. Three are included in this package.
#'
#' @param expt Experiment level phyloseq object.
#' @param taxrank Rank for which to make comparisons.
#' @param grps  Factor in sample data for which to make comparisons.
#' @param transformation Transformation function to use.
#' @param pc.filter Minimum percentage of total counts to include rank in result.
#' @param p.adjust.method Adjustment method for multiple comparisons.
#' @param pool.sd Logical, whether or not to pool standard deviations.
#'
#' @details transformation may be "none" or a user-supplied function name in quotation marks or any of the built-it transformations ("arc_sine", "log_arc_sine", or "sqrt_arc_sine"). The "sqrt_arc_sine" has generaally proven most effective.
#'
#' @return A data frame with taxa as rows and results in columns.
#'
#' @export
#'
#' @importFrom multcompView multcompLetters
#'
#' @seealso arc_sine, log_arc_sine, sqrt_arc_sine, check_var
#'
#' @examples
#'
#'
make_comparisons <- function(expt, taxrank = "Phylum",  grps = "Treatment", transformation = "none", pc.filter = 0.01, p.adjust.method ="BH",  pool.sd = FALSE) {

  # Prepare phyloseq objects
  expt.pc <- transform_sample_counts(expt, function(x) 100*(x/sum(x)))

  # Agglomerate to desired rank.
  expt.taxon <- tax_glom(expt, taxrank)
  expt.taxon.pc <- tax_glom(expt.pc, taxrank)

  # Remove ranks other than taxrank
  tax_table(expt.taxon) <- tax_table(expt.taxon)[ , taxrank]
  tax_table(expt.taxon.pc) <- tax_table(expt.taxon.pc)[ , taxrank]

  # Filter out taxa that are < 0.1% of the total sequences in expt.
  n <- sum(taxa_sums(expt)) * pc.filter
  expt.taxon <- prune_taxa(taxa_sums(expt.taxon)>=n, expt.taxon)
  expt.taxon.pc <- prune_taxa(taxa_names(expt.taxon), expt.taxon.pc)

  # Extract percentage otu table, make rownames the taxa names.
  otu.pc <- veganotu(expt.taxon.pc)
  colnames(otu.pc) <- tax_table(expt.taxon.pc)[colnames(otu.pc), taxrank]

  # Apply transformation
  if (!(transformation == "none")) {
    otu.pc.trans <- apply(otu.pc, MARGIN = 1, FUN = transformation)
  } else {
    otu.pc.trans <- t(otu.pc)
  }

  otu.pc <- t(otu.pc)

  # Extract sample data table, get vector of groups
  sam <- vegansam(expt.taxon.pc)
  grp <- sam[ , grps]
  grp <- factor(grp, ordered = FALSE)

  # Calculate mean and standard deviation for each row of data.
  otu.pc.means <- apply(otu.pc, 1, mean)
  otu.pc.sd <- apply(otu.pc, 1, sd)
  otu.pc.sum <- data.frame(otu.pc.means, otu.pc.sd)
  colnames(otu.pc.sum) <- c("mean", "sd")
  otu.pc.sum <- round(otu.pc.sum, 2)

  # Make global F tests
  # Prepare result data frame
  rslt <- as.data.frame(matrix(nrow=nrow(otu.pc), ncol=3))
  colnames(rslt) <- c(taxrank, "F", "Prob")
  rslt[ , 1] <- rownames(otu.pc)

  for (i in 1:nrow(otu.pc.trans)) {
    temp2 <- oneway.test(as.numeric(otu.pc.trans[ i, ]) ~ grp, var.equal = pool.sd)
    rslt[i,2] <- temp2$statistic
    rslt[i,3] <- temp2$p.value
  }

  rslt$sig <- lapply(rslt$Prob, asterix)
  rslt$F_value <- paste(round(rslt$F, 2), rslt$sig, sep="")
  rownames(rslt) <- rslt[ , 1]
  f.test.rslt <- rslt[ , -1]

  part1 <- merge_2_frames(otu.pc.sum, f.test.rslt)
  part1 <- part1[ , -c(3:5)]

  # Prepare comparison results matrix
  comp.sum <- matrix(NA, nrow=nrow(otu.pc), ncol=(length(levels(grp)) + 1))
  colnames(comp.sum) <- c(taxrank, levels(grp))

  # Fill comparison result matrix
  for (i in 1:nrow(otu.pc.trans)) {
    taxrank.name <- as.character(rownames(otu.pc)[i]) # gets taxon name

    # Calculate means and standard deviations by row
    df <- data.frame(otu.pc[i, ], grp)
    colnames(df) <- c("value", grps)
    avg <- tapply(df[ , 1], df[ , 2], mean) # calculates means
    std.dev <- tapply(df[ , 1], df[ , 2], sd) # calculates sd
    avg <- round(avg, 2)
    std.dev <- round(std.dev, 2)
    l <- paste(avg, "+/-", std.dev, sep = "")

    # Pairwise t-tests
    x <- pairwise.t.test(otu.pc.trans[ i, ], g = grp, p.adjust.method = p.adjust.method, pool.sd = pool.sd)

    # Get letter assignments for results; allow for NAs, NaNs in results
    x <- get_groups(x)$p.matrix
    x <- replace(x, is.na(x), 0)
    ltrs <- multcompView::multcompLetters(x)
    ltrs <- unname(ltrs$Letters)

    # Paste results together
    comp <- paste(l, ltrs, sep="")
    comp.sum[ i, ] <- c(taxrank.name, comp)
  }

  # Make first column rownames, then drop first column
  rownames(comp.sum) <- comp.sum[ , 1]
  comp.sum <- comp.sum[ , -1]

  comp.table <- merge_2_frames(part1, comp.sum)
  colnames(comp.table)[3] <- "F"

  return(list(comparison.table = comp.table, taxa.pc = otu.pc, taxa.pc.transformed = otu.pc.trans, groups = grp))

}
