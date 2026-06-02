#' Extract FastQC Adapter Content Table
#'
#' Parses a FastQC `.zip` file to locate the "Adapter Content" module and extracts 
#' the data as a clean data frame, preserving read positions as ordered categories.
#'
#' @param zip_path A character string specifying the path to the FastQC `.zip` file.
#'
#' @return A `data.frame` where the first column (`Position`) is an ordered factor,
#'   subsequent columns represent the cumulative percentages of different adapter types,
#'   and two trailing columns: `sample_id` (derived from the file name by stripping
#'   FastQC suffixes and read-direction tags) and `direction` (`"forward"`, `"reverse"`,
#'   or `NA` if the direction cannot be determined from the file name). Returns `NULL`
#'   if the section is missing or empty.
#' @export
extract_adapter_content <- function(zip_path) {
  
  file_name <- basename(zip_path)
  
  # 1. List files inside the zip to find fastqc_data.txt
  zip_files <- utils::unzip(zip_path, list = TRUE)$Name
  data_file <- grep("fastqc_data\\.txt$", zip_files, value = TRUE)
  
  if (length(data_file) == 0) {
    stop(paste("Error: fastqc_data.txt not found in", file_name))
  }
  
  # 2. Read the lines directly from the zipped file connection
  con <- unz(zip_path, data_file)
  lines <- readLines(con, warn = FALSE)
  close(con)
  
  # 3. Locate the boundaries of the Adapter Content module
  start_idx <- grep("^>>Adapter Content", lines)
  if (length(start_idx) == 0) {
    warning(paste("Section 'Adapter Content' not found in", file_name))
    return(NULL)
  }
  
  all_ends <- grep("^>>END_MODULE", lines)
  end_idx <- all_ends[all_ends > start_idx][1]
  
  section_lines <- lines[(start_idx + 1):(end_idx - 1)]
  
  # 4. Extract and clean the data
  header_line <- grep("^#Position", section_lines, value = TRUE)
  if (length(header_line) == 0) return(NULL)
  
  col_names <- unlist(strsplit(gsub("^#", "", header_line), "\t"))
  data_lines <- section_lines[!grepl("^#", section_lines)]
  
  if (length(data_lines) == 0) return(NULL)
  
  # 5. Reconstruct data lines into a clean data frame
  split_data <- strsplit(data_lines, "\t")
  data_matrix <- do.call(rbind, split_data)
  
  adapter_df <- as.data.frame(data_matrix, stringsAsFactors = FALSE)
  colnames(adapter_df) <- col_names
  
  # 6. Convert numerical adapter data columns
  if (ncol(adapter_df) > 1) {
    for (i in 2:ncol(adapter_df)) {
      adapter_df[[i]] <- as.numeric(adapter_df[[i]])
    }
  }
  
  # 7. CRITICAL: Convert Position into an Ordered Factor based on original file order
  # This prevents R/ggplot from sorting "100-104" before "2" alphabetically.
  adapter_df$Position <- factor(adapter_df$Position, 
                                levels = unique(adapter_df$Position), 
                                ordered = TRUE)
  
  # 8. Derive sample ID by stripping FastQC suffix and read-direction tags
  sample_id <- sub("_fastqc$", "", tools::file_path_sans_ext(file_name))
  sample_id <- sub("_R?[12]([^0-9].*)?$", "", sample_id)
  
  # 9. Derive read direction from file name
  direction <- if (grepl("_R?1([^0-9]|$)", file_name)) {
    "forward"
  } else if (grepl("_R?2([^0-9]|$)", file_name)) {
    "reverse"
  } else {
    NA_character_
  }
  
  adapter_df <- cbind(sample_id = sample_id, direction = direction, adapter_df)
  
  return(adapter_df)
}
