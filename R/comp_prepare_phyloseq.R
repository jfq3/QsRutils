#' Prepare Phyloseq
#'
#' Prepares a phyloseq object for making comparisons of relative abundances among treatments.
#'
#' @param expt Experiment level phyloseq object.
#' @param taxrank Taxonomic rank for which to make comparisons.
#' @param pc.filter Minimum percentage of total counts to include rank in result.
#'
#' @return A list of two modified experiemnt level phyloseq objects
#'
#' @details The otu_table in one of the returned objects has been transformed to percentages based on the original phyloseq object supplied. The taxa in both have been filtered to include only OTUs initially present at >= pc.filter times the original total counts. For both only taxrank is included in the tax_table.
#'
#' @export
#'
#' @examples
#'
#'
comp_prepare_phyloseq <- function(expt, taxrank = "Phylum", pc.filter = 0.01) {
  # Make copy with percentages instead of counts.
  expt.pc <- transform_sample_counts(expt, function(x) 100*(x/sum(x)))

  # Agglomerate to desired rank.
  expt.taxon <- tax_glom(expt, taxrank)
  expt.taxon.pc <- tax_glom(expt.pc, taxrank)

  # Remove ranks other than taxrank.
  tax_table(expt.taxon) <- tax_table(expt.taxon)[ , taxrank]
  tax_table(expt.taxon.pc) <- tax_table(expt.taxon.pc)[ , taxrank]

  # Filter out taxa that are < 0.1% of the total sequences in expt.
  n <- sum(taxa_sums(expt)) * pc.filter
  expt.taxon <- prune_taxa(taxa_sums(expt.taxon)>=n, expt.taxon)
  expt.taxon.pc <- prune_taxa(taxa_names(expt.taxon), expt.taxon.pc)

  return(list(expt.taxon = expt.taxon, expt.taxon.pc = expt.taxon.pc))
}
