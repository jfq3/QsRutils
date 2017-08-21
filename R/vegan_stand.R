#' Standardize a Phyloseq OTU Table
#'
#' Applies any vegan decostand standardization method to a phyloseq OTU table.
#'
#' @param physeq A phyloseq object containing at least an OTU table.
#' @param method A method from vegan's decostand function.
#' @param ... Other parameters passed to vegan's decostand function.
#'
#' @return Returns a phyloseq object with transformed OTU table.
#' @export
#'
#' @examples
#'
#' @importFrom phyloseq taxa_are_rows
#' @importClassesFrom phyloseq otu_table
#' @importFrom vegan decostand
#'
vegan_stand <-
function(physeq, method="hellinger", ...) {
  test <- taxa_are_rows(physeq)
  OTU <- otu_table(physeq)
  OTU <- as(OTU, "matrix")
  if (test) {
    OTU <- t(OTU)
  }
  OTU <- decostand(OTU, method, ...)
  if (test) {
    OTU <- t(OTU)
  }
  otu_table(physeq) <- otu_table(OTU, taxa_are_rows=test)
  return(physeq)
}
