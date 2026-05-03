
#' Subset physeq by refseq lengths
#'
#' @param p An experiment level phyloseq object with reference sequences
#' @param min_len The minimum reference sequence length to keep
#' @param max_len The maximum references sequences length to keep
#'
#' @returns A phyloseq object filtered to have reference sequences within the length range specified.
#' @details Sometimes, due to sequencing errors, reference sequences have lengths less than and/or greater than the expected amplicon length. This function offers an easy way of removing such extraneous reference sequences from an experiment level phyloseq object. The default values for min_len and max_len are fro the V4 region of the 16S rRNA gene.
#' @export
#'
#' @examples
#' {
#' data("expt")
#' expt@refseq@ranges@width |> 
#'   summary()
#' expt_filt <- subset_by_refseq_lengths(p = expt,
#'                                       min_len = 400,
#'                                       max_len = 420)
#' expt_filt@refseq@ranges@width |> 
#'   summary()
#' }
#' 
subset_by_refseq_lengths <- function(p, min_len = 252, max_len = 255) {
  # Get reference sequences
  refseqs <- phyloseq::refseq(p)
  if (is.null(refseqs)) {
    stop("No reference sequences found in the phyloseq object 'p'")
  }
  
  # Compute lengths: prefer Biostrings::width for XStringSet objects
    seq_lengths <- Biostrings::width(refseqs)

  # Build a tibble with explicit column names to avoid name collisions (don't use 'length')
  df <- tibble::tibble(
    taxon = names(refseqs),
    seq_length = as.integer(seq_lengths)
  )
  
  taxa2keep <- df |>
    dplyr::filter(seq_length >= as.integer(min_len) & seq_length <= as.integer(max_len)) |>
    dplyr::pull(taxon)
  
  p_out <- phyloseq::prune_taxa(taxa2keep, p)
  p_out
}

