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
#'
#'
comp_assemble <- function(part1, part2, part3) {
  comp <- merge_2_frames(part1, part2)
  comp <- comp[ , -c(3:5)]
  comp <- merge_2_frames(comp, part3)
  return(comp)
}
