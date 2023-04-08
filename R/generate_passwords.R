#' Genererate a Password
#'
#' Generates a random character string of specified length.
#'
#' @param n Numberof characters in password.
#' @param type c("alpha_numeric", "anything_else")
#'
#' @return A character string.
#' @export
#' @details If type equals "alpha_numeric" (the default), only alpha-numeric characters are used to generate the password. If type does not equal "alpha_numeric" then at least one non-alpha-numeric symbol will be included in the password. In either case, the alpha characters used are both upper and lower case.
#'
#' @examples
#' generate_password(8)
#'
generate_password <- function(n, type = "alpha_numeric") {
  set.seed(as.integer(Sys.time()))
  symb <- c("!", "@", "#", "$", "%", "^", "&", "?")
  let.nos <- c(letters, toupper(letters), as.character(0:9))
  all.symbols <- c(let.nos, symb, symb)
  pw <- ""
  if (type == "alpha_numeric") {
    pw <- sample(let.nos, n, replace=TRUE)
  } else {
    while (!(any(symb %in% pw))) {
      pw <- sample(all.symbols, n, replace=TRUE)
    }
  }
  pass <- NULL
  for (i in 1:n) {
    pass <- paste(pass, pw[i], sep="")
  }
  return(pass)
}
