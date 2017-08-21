#' Clear Warnings
#'
#' Clears all warning messages from the base environment.
#'
#' @export
#'
#' @details  Sometimes when working in the console R retains a list of warnings such that they keep being reported after the function call which originated them. This function removes them so that they are not a nuisance
#'
#' @examples
#'
#' clear_warnings()
#'
clear_warnings <- function() {
  assign("last.warning", NULL, envir = baseenv())
}
