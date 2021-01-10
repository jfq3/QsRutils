#' Root Tree in phyloseq Object
#'
#' Roots an unrooted tree in a phyloseq object
#' 
#' @param phylo A phyloseq object containing an unrooted tree
#' @details The tree is rooted by the longest terminal branch.
#' @return The same phyloseq object with a rooted tree
#' @export
#' @importFrom  ape Ntip
#' @importFrom magrittr %>%
#' @importFrom data.table data.table
#' @examples
#' \dontrun{
#' expt.rooted <- root_phyloseq_tree(expt.unrooted)
#' }
root_phyloseq_tree <- function(phylo) {
  tree.unrooted <- phy_tree(phylo)
  # tablify parts of tree that we need.
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  new.tree <- ape::root(tree.unrooted, outgroup=new.outgroup, resolve.root=TRUE)
  phy_tree(phylo) <- new.tree
  return(phylo)
}
