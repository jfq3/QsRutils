#' @title Remove elements from a vector
#'
#' @description This infix operator removes all elements of `y` from `x`.
#'
#' @param x A vector.
#' @param y A vector of elements to remove from `x`.
#'
#' @return A vector containing elements from `x` that are not present in `y`.
#'
#' @examples
#' c(1, 2, 3, 4, 5) %w/o% c(2, 4)
#' letters[1:5] %w/o% c("a", "c")
#'
#' @export
`%w/o%` <- function(x, y) {
  x[!x %in% y]
}

