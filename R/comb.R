#' comb
#'
#' Calculates the number of combinations of n things drawn r at a time.
#'
#' @param n The total number of items.
#' @param r The number of items to be drawn.
#' @param repetition A logical, whether or not repetitions are allowed. FALSE by default.
#' 
#' @return An integer giving the number of ways a set of r items can be drawn from a set of n items.
#' @export
#'
#' @examples
#' comb(5, 3)
#' comb(5, 3, repetition = TRUE)
#'
comb = function(n, r, repetition = FALSE) {
  if (repetition != TRUE) {
    return(factorial(n) / (factorial(r) * factorial(n-r)))
  } else {
    return(factorial(r+n-1)/(factorial(r)*factorial(n-1)))
  }
}
