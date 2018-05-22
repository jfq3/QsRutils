#' Permutations
#'
#' Retuns the number of permutaions of n things taken r at a time.
#'
#' @param n Total number of items.
#' @param r Number of items drawn.
#' @param repetition A logical, whether or not repetitions are allowed. FALSE by default.
#'
#' @return An integer giving how many ways m things can be drawn n at a time.
#' @export
#'
#' @examples
#' perm(10, 5)
#' perm(10, 5, repetition = TRUE)
#'
perm =   function(n, r, repetition = FALSE) {
  if (repetition != TRUE) {
    return(factorial(n) / factorial(n-r))
  } else {
    return(n^r)
  }
}
