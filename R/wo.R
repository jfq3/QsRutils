#' @title Remove elements from a vector
#' @description Removes elements present in `y` from `x`.
#' @param x A vector.
#' @param y A vector of elements to remove from `x`.
#' @return A vector containing elements of `x` that are not in `y`.
#' @examples
#' c(1, 2, 3, 4, 5) %wo% c(2, 4)
#' letters[1:5] %wo% c("b", "d")
#' @export
`%wo%` <- function(x, y) {
  x[!x %in% y]
}