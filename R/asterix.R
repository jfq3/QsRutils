#' Indicate Significance with Stars
#'
#' @param prob p value
#'
#' @return Character vector of asterisks indicating significance level.
#' @export
#'
#' @details Returns '***; for p < 0.001, '**' for p < 0.01, '*' for p < 0.05.
#' @examples
#' asterix(0.039)
#'
asterix <- function(prob) {
  if (is.na(prob)) {
    sig <- c("")
  }
  else if(prob<0.001) {
    sig <- c("***")
  }
  else if(prob<0.01) {
    sig <- c("**")
  }
  else if(prob<0.05) {
    sig <- c("*")
  }
  else {
    sig<- c("")
  }
  return(sig)
}
