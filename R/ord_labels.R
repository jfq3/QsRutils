#'Make Ordination Axis Labels
#'
#'Makes ordination axis labels that include, if appropriate, the % of the total
#'variance explained by each axis.
#'
#'@usage ord_labels(ord)
#'
#'@param ord A vegan ordination object.
#'
#'@return A character vector, each element of which can be used to label the
#'  corresponding axis of an ordination plot.
#'
#'@details 
#'  If there are no eigenvalues in ord, or if any eigenvalues are less than 0,
#'  each element of the vector returned has the form "DIMn" where n is the axis
#'  number. Otherwise, each element of the vector returned has the form Axisn
#'  xx.x % where Axis is taken from the vector of eigenvalues in ord if they
#'  are named or simply "DIM" if they are not, n is the number of the axis, and
#'  xx.x is the % of total variance explained by the axis.
#'  
#'  For this function to work correctly, ord should be created in one of
#'  the following ways: 
#'  1. As an unconstrained ordination using vegan::rda. In
#'  this case the labels have the form PCAn xx.x % 
#'  2. As a PCoA made with stats::cmdscale. In this case the labels have the form DIMn xx.x %.
#'  3. As a CA made with vegan::ca. In this case the labels have the form CAn xx.x %.
#'  
#'@md
#'@seealso [rda_labels()]
#'@export
#'
#'@examples
#'{
#' # For PCA using rda:
#' data("dune", package = "vegan")
#' dune_hel <- vegan::decostand(dune, method = "hellinger")
#' pca <- vegan::rda(dune_hel)
#' ord_labels(pca)[1:2]
#' 
#' # For PCoA with negative eigenvalues
#' d <- vegan::vegdist(dune)
#' pcoa <- stats::cmdscale(d, k = nrow(dune)-1, eig = TRUE, add = FALSE)
#' ord_labels(pcoa)[1:2]
#' 
#' # For PCoA without negative eigenvalues
#' pcoa <- stats::cmdscale(d, k = nrow(dune)-1, eig = TRUE, add = TRUE)
#' ord_labels(pcoa)[1:2]
#' 
#' # For correspondence analysis
#' ca_ord <- vegan::ca(dune)
#' ord_labels(ca_ord)
#'}

ord_labels <-
  function(ord){
    ev <- vegan::eigenvals(ord)
    tol <- (1e-07)*ev[1]
    ord.labels <- rep("", length(ev))
    if (any(ev < -tol)) {
      for ( i in 1:length(ev)) {
        ord.labels[i] <- paste("DIM", i, sep = "")
      }
    } else {
      ev.pc <- round(100*(ev/sum(ev)), 2)
      axis.names <- names(ev)
      if (is.null(axis.names)) {
        for ( i in 1:length(ev.pc)) {
          ord.labels[i] <- paste("DIM", i, " ", sprintf(ev.pc[i], fmt = '%#.1f'), "%", sep="")
        }
      } else {
        for (i in 1:length(ev.pc)){
          ord.labels[i] <- paste(axis.names[i], " ", ev.pc[i],"%", sep="")
        } 
      }
    }
    return(ord.labels)
  }
