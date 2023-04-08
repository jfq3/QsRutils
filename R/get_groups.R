#' get_groups
#'
#' Assign treatment groups based on pairwise t-tests.
#'
#' @param ptt.rslt Result from the stats function pairwise.t.test.
#' @param alpha Confidence level.
#' @param rm.subset A logical; remove group subsets if true.
#'
#' @return A list consisting of groups of treatment groups that are not significantly differnet and a matrix of p values.
#'
#' @details This function aids in making letter assignments as to which treatments are significantly different. Also returns a square matrix of alpha values for all pairwise differences. This square matrix can serve as input to the multcompLetters function of the multcompView package which provides letter assignments.
#' If rm.subset is FALSE, then groups such as {A,B} and {A, B, C} may be reported. This is redundant in the sense the {A, B} is a subset of {A, B, C}. In this case if rm.subset is FALSE, the grop {A, B} is not reported.
#'
#' @seealso make_letter_assignments
#' @export
#'
#' @examples
#' \dontrun{
#' get_groups(ptt.rslt, alpha = 0.05, rm.subset = FALSE)
#' }
#'
get_groups <- function(ptt.rslt, alpha = 0.05, rm.subset = FALSE) {

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

  # For each column, get list of treatments not significantly different.
  grp <- list()
  trts <- colnames(my.mat)
  for (i in 1:ncol(my.mat)) {
    grp[[i]] <- c(trts[i], names(which(my.mat[ , i] > alpha)))
  }

  # Remove groups that are sub-sets of other groups
  k <- 0
  del <- vector()
  for (i in 1:(length(grp)-1)) {
    for ( j in (i+1):length(grp)) {
      if (!rm.subset) {
        if (setequal(grp[[i]], grp[[j]])) {
          k <- k+1
          del[k] <- j
        }
      }
      else {
        if (all(is.element(grp[[i]], grp[[j]]))) {
          k <- k+1
          del[k] <- i
        }
        else if (all(is.element(grp[[j]], grp[[i]]))) {
          k <- k+1
          del[k] <- j
        }
      }

    }
  }

  del <- unique(del)
  del <- del[order(del, decreasing = TRUE)]

  if (length(del) >= 1) {
	  for (i in 1:length(del)) {
		grp[[del[i]]] <- NULL
	  }
  }

  return(list(groups = grp, p.matrix = my.mat))

}
