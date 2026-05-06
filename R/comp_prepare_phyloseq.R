#' Prepare Phyloseq
#'
#' Prepares a phyloseq object for making comparisons of relative abundances
#' among treatments.
#'
#' @param expt Experiment level phyloseq object.
#' @param taxrank Taxonomic rank for which to make comparisons.
#' @param pc.filter Minimum percentage of total counts to include rank in
#'   result.
#'
#' @return A list of two modified experiment level phyloseq objects.
#'
#' @details The otu_table in one of the returned objects has been transformed to
#'   percentages based on the original phyloseq object supplied. The taxa in
#'   both have been filtered to include only OTUs initially present at >=
#'   pc.filter times the original total counts. For both only, taxrank is
#'   included in the tax_table.  
#'   See also the vignette "Compare Relative Abundances Among Treatments."
#'   
#'      
#' @export
#' @importFrom phyloseq prune_taxa
#' @importFrom phyloseq tax_glom
#' @importFrom phyloseq taxa_names
#' @importFrom phyloseq taxa_sums
#' @importFrom phyloseq tax_table
#' @importFrom phyloseq transform_sample_counts
#' @examples
#' {
#' data(its.root)
#' temp1 <- comp_prepare_phyloseq(its.root)
#' temp1
#' }
#' 
comp_prepare_phyloseq <- function(expt, taxrank = "Phylum", pc.filter = 0.01) {
  # Make copy with percentages instead of counts.
  expt.pc <- phyloseq::transform_sample_counts(expt, function(x) 100*(x/sum(x)))

  # Agglomerate to desired rank.
  expt.taxon <- phyloseq::tax_glom(expt, taxrank)
  expt.taxon.pc <- phyloseq::tax_glom(expt.pc, taxrank)

  # Remove ranks other than taxrank.
  phyloseq::tax_table(expt.taxon) <- phyloseq::tax_table(expt.taxon)[ , taxrank]
  phyloseq::tax_table(expt.taxon.pc) <- phyloseq::tax_table(expt.taxon.pc)[ , taxrank]

  # Filter out taxa that are < pc.filter of the total sequences in expt.
  n <- sum(taxa_sums(expt)) * pc.filter
  expt.taxon <- phyloseq::prune_taxa(phyloseq::taxa_sums(expt.taxon)>=n, expt.taxon)
  expt.taxon.pc <- phyloseq::prune_taxa(phyloseq::taxa_names(expt.taxon), expt.taxon.pc)

  return(list(expt.taxon = expt.taxon, expt.taxon.pc = expt.taxon.pc))
}
