#' Merge Two Data Frames
#'
#' Merge two data frames by their row names.
#'
#' @param one A data frame.
#' @param two A second data frame.
#'
#' @return A merged data frame.
#' @export
#' @details Merges data frames by common row names. This function differs from  merge.data.frames in that the merged data frame returned has row names and not a new column of the row names.
#' @examples
#' {
#' common_rows <- paste0("ID_", 1:5)
#'
#' df1 <- data.frame(
#'   Value_A = runif(5), 
#'   Category = sample(c("X", "Y"), 5, replace = TRUE),
#'   row.names = common_rows
#' )
#' 
#' df2 <- data.frame(
#'   Value_B = rnorm(5),
#'   Flag = sample(c(TRUE, FALSE), 5, replace = TRUE),
#'   row.names = common_rows
#' )
#' 
#' merge_2_frames(df1, df2)
#' }
#'
merge_2_frames <- function(one, two) {
  rslt <- merge.data.frame(one, two, by = 0)
  rownames(rslt) <- rslt[ , 1]
  rslt <- rslt[ , -1]
  return(rslt)
}
