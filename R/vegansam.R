#' Extract Sample Data Table
#'
#' Extracts a sample data table from a phyloseq object.
#'
#' @param physeq A phyloseq object containing sample_data.
#'
#' @return A data frame with samples in rows and factors and/or variables in columns.
#' @export
#'
#' @examples
#' \dontrun{
#' vegansam(physeq)
#' }
#'
#' @importFrom phyloseq sample_data
#' @importFrom phyloseq nsamples
#' @importFrom phyloseq sample_variables
#' @importFrom phyloseq sample_names
#'
vegansam <- function(physeq) {
  sam <- sample_data(physeq)
  c <- length(sample_variables(physeq))
  r <- nsamples(physeq)
  sam2 <- as.data.frame(matrix(data = NA, nrow = r, ncol = c))
  for (i in 1:c) {
    sam2[ , i] <- sam@.Data[[i]]
  }
  rownames(sam2) <- sample_names(physeq)
  colnames(sam2) <- sample_variables(physeq)
  return(sam2)
}
