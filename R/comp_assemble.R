#' Assemble Comparison Parts
#'
#' Assembles Comparison Data Frame
#'
#' @param part1 Result from comp_means_sd
#' @param part2 Result from comp_make_f_tests
#' @param part3 Result from comp_comparisons
#'
#' @return A summary data frame of differential abundances by taxon and treatment.
#' @export
#'
#' @examples
#' {
#' data("its.root")
#' temp1 <- comp_prepare_phyloseq(its.root)
#' temp2 <- comp_prepare_otu_table(temp1$expt.taxon.pc,
#'                                 grps = "Label",
#'                                 transformation = "sqrt_arc_sine")
#' temp3 <- comp_means_sd(temp2$otu.pc)
#' temp4 <- comp_make_f_tests(temp2$otu.pc.trans,
#'                            grps = temp2$groups,
#'                            var.equal = TRUE)
#' temp5 <- comp_comparisons(otu.pc = temp2$otu.pc,
#'                           otu.pc.trans = temp2$otu.pc.trans,
#'                           grps = temp2$groups,
#'                           p.adjust.method = "BH",
#'                           pool.sd = TRUE)
#' comp_assemble(temp3, temp4, temp5) |> 
#'   dplyr::arrange(desc(mean))
#' }
#'
comp_assemble <- function(part1, part2, part3) {
  comp <- merge_2_frames(part1, part2)
  comp <- comp[ , -c(3:5)]
  comp <- merge_2_frames(comp, part3)
  return(comp)
}
