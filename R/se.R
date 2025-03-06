#' Calculate the standard error
#' 
#' Calculate the standard error of a numeric vector allowing for NA values.
#' 
#' @usage se(x)
#' @param x A numeric vector
#' @return The standard error
#' @export
#' @details Curiously there is no R function for calculating the standard error. This function allows for NA values.
#' @author John Quensen
#' 
#' se <- function(x) {
#'  se <-  sd(x) / sqrt(sum(!is.na(x)))
#'  return(se)
#' }
#' 
#' @examples
#' #' 
#' x <- c(2, 6, 4, 5, 2, NA)
#' se = sd(x) / sqrt(sum(!is.na(x)))
#' se
#' 
#' 