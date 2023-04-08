#' Extract Vegan OTU Table
#'
#' Extracts a vegan compatible OTU table from a phyloseq object.
#'
#' @param physeq A phyloseq object contaning at least an OTU table.
#'
#' @return A matrix with samples in rows and OTUs in columns.
#' @export
#'
#' @example
#' \dontrun{
#' veganotu(physeq)
#' }
#'
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq taxa_are_rows
#'
veganotu <-
function(physeq) {
    OTU <- otu_table(physeq)
    if (taxa_are_rows(OTU)) {
        OTU <- t(OTU)
    }
    OTU <- as(OTU, "matrix")
    return(OTU)
}
