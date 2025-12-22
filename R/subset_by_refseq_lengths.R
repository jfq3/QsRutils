
#' Subset physeq by refseq lengths
#'
#' @param p An experiment level phyloseq object with reference sequences
#' @param min_len The minimum reference sequence length to keep
#' @param max_len The maximum references sequences lengtht to keep
#'
#' @returns A phyloseq object filtered to have reference sequences within the length range specified.
#' @details Sometimes, due to sequencing errors, reference sequences have lengths less than and/or greater than the expected lenght of the expected amplicon length. This function offeres an easy way of removing such extraneous reference sequences from an experiment level phyloseq object. Teh default values of min_len and amx_len are fro the V4 region of the 16S rRNA gene.
#' @export
#'
#' @examples
#' \dontrun {
#' p <- subset_by_refseq_lengths(p)
#' }
subset_by_refseq_lengths <- function(p, min_len = 252, max_len = 255) {
  refseqs_data <- refseq(p) 
  taxa2keep <- nchar(as.character(refseqs_data)) |> 
    as.data.frame() |> 
    rownames_to_column(var = "taxon") |> 
    rename(length = `nchar(as.character(refseqs_data))`) |> 
    filter(length >= min_len & length <= max_len) |> 
    pull(taxon)
  p <- prune_taxa(taxa2keep, p)
  return(p)
}

