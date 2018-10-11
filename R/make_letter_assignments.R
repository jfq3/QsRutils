#' Make Letter Assignments
#'
#' Makes letter assignments for treatment groups that are not significantly different.
#'
#' @param ptt.rslt Output from the pairwise.t.test function.
#' @param significance Alpha level to be declared a significant difference.
#'
#' @return Lists of letter assignments.
#'
#' @details Letter assignments are made using Piepho's algorithm.
#'
#' @references Piepho, H. P. 2004. An algorithm for a letter-based representation of all-pairwise comparisons. Journal of Computational and Graphical Statistics **13**:456-466.
#'
#' @export
#'
#' @examples
#'
#' @importFrom multcompView multcompLetters
#'
make_letter_assignments <- function(ptt.rslt, significance=0.05) {

  # Get a matrix of the alpha values
  temp <- ptt.rslt$p.value

  # Make a square matrix to populate with the alpha values.
  n <- nrow(temp)
  mat.names <- c(colnames(temp), rownames(temp)[n])
  my.mat <- matrix(data = NA, nrow = n+1, ncol = n+1)
  colnames(my.mat) <- mat.names
  rownames(my.mat) <- mat.names

  # Add diagonal.
  for (i in 1:nrow(my.mat)) {
    my.mat[ i, i] <- 0
  }

  # Get vector of p.values
  stat <- na.exclude(as.vector(ptt.rslt$p.value))

  # Add other cells to square matrix.
  k=1
  for (j in 1:(nrow(my.mat)-1)) {
    for (i in ((j+1):nrow(my.mat))) {
      my.mat[i,j] <-  my.mat[j,i] <- stat[k]
      k=k+1
    }
  }
  letter.assignments <- multcompView::multcompLetters(my.mat, threshold=significance)
  return(letter.assignments)
}
