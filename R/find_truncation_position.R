#' Find Read Truncation Position Based on FastQC Quality and Read Direction
#'
#' Parses a FastQC `.zip` file for a specific read direction (forward or reverse), 
#' looks for the "Per base sequence quality" section, and identifies the first 
#' position where the Lower Quartile dips to or below a specified Q-score threshold.
#'
#' @param zip_path A character string specifying the path to the FastQC `.zip` file.
#' @param read_direction A character string, either `"forward"` or `"reverse"`. 
#'   Matches common patterns like `_R1_`/`_1` for forward and `_R2_`/`_2` for reverse.
#' @param target_q A numeric value specifying the Lower Quartile Q-score threshold. Default is 20.
#'
#' @return A one-row data frame with columns `sample_name` (file name without `_fastqc.zip`),
#'   `direction` (`"forward"` or `"reverse"`), and `position` (quality result string).
#'   Returns `NULL` if the file does not match the requested read direction.
#' @details
#'   Results for all files in a path can be aggregated with \code{purrr::list_rbind()}:
#'   \preformatted{
#'   all_zips <- list.files(path = "path_to_zip_files", pattern = "\\.zip$", full.names = TRUE)
#'   purrr::map(all_zips, find_truncation_position, "forward") |> purrr::list_rbind()
#'   purrr::map(all_zips, find_truncation_position, "reverse") |> purrr::list_rbind()
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' find_truncation_position("sample1_L001_R1_001_fastqc.zip", read_direction = "forward")
#' }
find_truncation_position <- function(zip_path, read_direction = c("forward", "reverse"), target_q = 20) {
  
  # Match and validate the read direction argument
  read_direction <- match.arg(read_direction)
  file_name <- basename(zip_path)
  
  # Define regex patterns for forward vs reverse reads
  is_forward <- grepl("(_R1_|_1(_|\\.fastqc))", file_name, ignore.case = TRUE)
  is_reverse <- grepl("(_R2_|_2(_|\\.fastqc))", file_name, ignore.case = TRUE)
  
  # Skip file if it doesn't match the user's selection
  if (read_direction == "forward" && !is_forward) return(NULL)
  if (read_direction == "reverse" && !is_reverse) return(NULL)
  
  # 1. List files inside the zip to find fastqc_data.txt
  zip_files <- utils::unzip(zip_path, list = TRUE)$Name
  data_file <- grep("fastqc_data\\.txt$", zip_files, value = TRUE)
  
  if (length(data_file) == 0) {
    stop(paste("Error: fastqc_data.txt not found in", file_name))
  }
  
  # 2. Read the specific file connections directly from the zip
  con <- unz(zip_path, data_file)
  lines <- readLines(con, warn = FALSE)
  close(con)
  
  # 3. Locate the target section boundaries
  start_idx <- grep("^>>Per base sequence quality", lines)
  if (length(start_idx) == 0) {
    return("Section 'Per base sequence quality' not found.")
  }
  
  all_ends <- grep("^>>END_MODULE", lines)
  end_idx <- all_ends[all_ends > start_idx][1]
  
  section_lines <- lines[(start_idx + 1):(end_idx - 1)]
  data_lines <- section_lines[!grepl("^#", section_lines)]
  
  # 4. Process line by line to locate the quality drop
  result <- paste0("Excellent quality: Lower Quartile never dips to Q", target_q)
  for (line in data_lines) {
    parts <- unlist(strsplit(line, "\t"))
    if (length(parts) < 5) next
    
    position <- parts[1]
    lower_quartile <- as.numeric(parts[4]) # 4th column is Lower Quartile
    
    if (!is.na(lower_quartile) && lower_quartile <= target_q) {
      result <- paste0("Dips to Q", lower_quartile, " at position ", position)
      break
    }
  }
  
  data.frame(
    sample_name = sub("_fastqc\\.zip$", "", file_name),
    direction   = read_direction,
    position    = result,
    stringsAsFactors = FALSE
  )
}
