
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
#' \dontrun{
#' p <- subset_by_refseq_lengths(p)
#' }
#' 
subset_by_refseq_lengths <- function(p, min_len = 252, max_len = 255) {
  `nchar(as.character(refseqs_data))` <- taxon <- NULL
  refseqs_data <- phyloseq::refseq(p) 
  taxa2keep <- base::nchar(base::as.character(refseqs_data)) |> 
    base::as.data.frame() |> 
    tibble::rownames_to_column(var = "taxon") |> 
    dplyr::rename(length = `nchar(as.character(refseqs_data))`) |> 
    dplyr::filter(length >= min_len & length <= max_len) |> 
    dplyr::pull(taxon)
  p <- phyloseq::prune_taxa(taxa2keep, p)
  return(p)
}

