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
#' {
#' data("expt")
#' print("Before standardization, first 5 rows and columns:")
#' print(veganotu(expt)[1:5, 1:5])
#' expt_mod <- vegan_stand(expt, method = "hellinger")
#' print("After standardization, first 5 rows and columns:")
#' print(veganotu(expt_mod)[1:5, 1:5])
#' }
#'
#' @importFrom phyloseq taxa_are_rows
#' @importClassesFrom phyloseq otu_table
#' @importFrom vegan decostand
#'
vegan_stand <-
function(physeq, method="hellinger", ...) {
  test <- phyloseq::taxa_are_rows(physeq)
  OTU <- phyloseq::otu_table(physeq)
  OTU <- as(OTU, "matrix")
  if (test) {
    OTU <- t(OTU)
  }
  OTU <- vegan::decostand(OTU, method, ...)
  if (test) {
    OTU <- t(OTU)
  }
  phyloseq::otu_table(physeq) <- phyloseq::otu_table(OTU, taxa_are_rows=test)
  return(physeq)
}
