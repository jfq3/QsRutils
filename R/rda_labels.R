#' Make RDA Axis Labels
#'
#' Makes RDA axis labels that include the % total variance explained by each axis.
#'
#' @param rda Object that contains CCA result from vegan's rda function.
#'
#' @return A character vector, each element of which can be used to label the corresponding axis of an RDA plot.
#'
#' @details Each element of the vector returned has the form "RDAn xx.x%" where n is the number of the RDA axis and xx.x is the % of total variance explained by the axis.
#' @export
#'
#' @examples
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
