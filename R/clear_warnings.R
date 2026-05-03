#' Clear Warnings
#'
#' Clears all warning messages from the base environment.
#'
#' @export
#'
#' @return No return value, called for side effects. Clears the stored warning
#'   list so that stale warnings are no longer reported in the console.
#'
#' @details Sometimes when working in the console R retains a list of warnings
#'   such that they keep being reported after the function call which originated
#'   them. This function removes them so that they are not a nuisance.
#'
#' @examples
#' {
#' as.numeric(c("1", "2", "apples"))
#' summary(warnings())
#' clear_warnings()
#' summary(warnings())
#' }
clear_warnings <- function() {
  tryCatch(
    assign("last.warning", NULL, envir = baseenv()),
    error = function(e) invisible(NULL)
  )
}
