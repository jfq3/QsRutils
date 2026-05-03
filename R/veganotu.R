#' Extract Vegan OTU Table
#'
#' Extracts a vegan compatible OTU table from a phyloseq object.
#'
#' @param physeq A phyloseq object contaning at least an OTU table.
#'
#' @return A matrix with samples in rows and OTUs in columns.
#' @export
#'
#' @examples
#' {
#' data("expt")
#' # Show only first 5 columns and rows:
#' veganotu(physeq = expt)[1:5, 1:5]
#' }
#'
#' @importFrom phyloseq otu_table
#' @importFrom phyloseq taxa_are_rows
#'
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
