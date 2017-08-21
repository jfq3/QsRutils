#' comb
#'
#' Calculates the number of combinations of n things drawn x at a time.
#'
#' @param n The total number of items.
#' @param x The number of items to be drawn.
#'
#' @return An integer giving the number of ways a set of x items can be drawn from a set of n items.
#' @export
#'
#' @examples
#'
#' comb(10, 2)
#'
comb = function(n, x) {
  return(factorial(n) / (factorial(x) * factorial(n-x)))
}
