#' Make RDA Axis Labels
#'
#' Makes RDA axis labels that include the % total variance explained by each RDA
#' axis, that is by the constrained portion of the analysis.
#'
#' @param rda A constrained ordination object made with vegan::rda() or
#'   vegan::cca().
#'
#' @return A character vector, each element of which can be used to label the
#'   corresponding axis of an RDA plot.
#'
#' @details Each element of the vector returned has the form "RDAn xx.x%" where
#'   n is the number of the RDA axis and xx.x is the % of total variance
#'   explained by the axis. The percent of total variance is for the constrained
#'   portion only.
#'
#' @export
#' @seealso [ord_labels()]
#' @examples
#' {
#' # Using vegan::rda()
#' data(dune, package = "vegan")
#' data(dune.env, package = "vegan")
#' dune.Manure <- vegan::rda(dune ~ Manure, dune.env)
#' rda_labels(dune.Manure)[1:2]
#' 
#' # Using vegan::cca()
#' data(varespec, package = "vegan")
#' data(varechem, package = "vegan")
#' vare.cca <- vegan::cca(varespec ~ Al + P*(K + Baresoil), data=varechem)
#' rda_labels(vare.cca)[1:2]
#' }
#' 
rda_labels <-
function(rda){
  ev <- rda$CCA$eig
  ev.pc <- round(100*(ev/sum(ev)), 2)
  rda.labels <- rep("a", length(ev.pc))
  for (i in 1:length(ev.pc)){
    rda.labels[i] <- paste("RDA",i," ", ev.pc[i],"%",sep="")
  }
return(rda.labels)
}
