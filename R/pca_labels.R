#' Make PCA Axis Labels
#'
#' Makes PCA axis labels that include the % total variance explained by each axis.
#'
#' @param pca Object containig the results of vegan's rda function.
#'
#' @return A character vector, each element of which can be used to label the corresponding axis of a PCA plot.
#'
#' @details Each element of the vector returned has the form "PCAn xx.x%" where n is the number of the PCA axis and xx.x is the % of total variance explained by the axis.
#'
#' @export
#'
#' @examples
#'
pca_labels <-
function(pca){
  ev <- pca$CA$eig
  ev.pc <- round(100*(ev/sum(ev)), 2)
  pca.labels <- rep("a", length(ev.pc))
  for (i in 1:length(ev.pc)){
    pca.labels[i] <- paste("PC",i," ", ev.pc[i],"%",sep="")
  }
return(pca.labels)
}
