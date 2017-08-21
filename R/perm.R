#' Permutations
#'
#' Retuns the number of permutaions of n things taken x at a time.
#'
#' @param n Total number of items.
#' @param x Number of items drawn.
#'
#' @return An integer giving how many ways n things can be drawn x at a time.
#' @export
#'
#' @examples
#' perm(5, 10)
#'
perm = function(n, x) {
  return(factorial(n) / factorial(n-x))
}
