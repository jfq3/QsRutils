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
#'
merge_2_frames <- function(one, two) {
  rslt <- merge.data.frame(one, two, by = 0)
  rownames(rslt) <- rslt[ , 1]
  rslt <- rslt[ , -1]
  return(rslt)
}
