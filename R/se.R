#' Standard error
#' 
#' Calculates the standard error of a numeric vector
#' 
#'
#' @param x A numeric vector
#'
#' @returns The standard error of the numeric vector
#' @details
#' NA values are ignored.
#' 
#' @export
#' 
#'
#' @examples
#' x <- c(1,2,3,4,5)
#' se(x)
#' 
se <- function(x) {
  se <-  sd(x) / sqrt(sum(!is.na(x)))
  return(se)
 }
