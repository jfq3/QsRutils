#' Calculate Good's Coverage
#'
#' Calculates Good's coverage from a community data matrix with samples as rows and OTUs as columns.
#'
#' @param com a vegan compatible community data matrix.
#'
#' @return A table with the headings number of singletons, number of sequences, and Good's coverage for each sample in rows.
#' @export
#'
#' @references Good, I. J. 1953. The Population Frequencies of Species and the Estimation of Population Parameters. Biometrika 40:237-264.
#'
#' @examples
#' \dontrun{
#' goods(species_matrix)
#' }
#'
goods <-
function(com){
  no.seqs <- rowSums(com)
  sing <- com==1
  no.sing <- apply(sing, 1, sum)
  goods <- 100*(1-no.sing/no.seqs)
  goods.sum <- cbind(no.sing, no.seqs, goods)
  goods.sum <- as.data.frame(goods.sum)
  return(goods.sum)
}
