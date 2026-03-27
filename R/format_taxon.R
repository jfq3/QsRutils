#' Format a taxon
#' 
#' Formats a taxon so that it will be properly italicized in a ggpplot
#'
#' @param x A string representing a phylum, class, etc.
#'
#' @returns The string with asterisks properly inserted.
#' @details If a taxon begins with an upper case letter followed by lower case letters and does not contain an underscore, it is wrapped in asterisks. If it begins with an upper case letter followed by lower case letters and contains an underscore, the portion before the underscore is wrapped in asterisks, the underscore removed and any letters following the underscore left alone. If a taxon contains all upper case letters or digits it is not a proper taxon and is left alone.
#' @export
#'
#' @examples
#' format_taxon("UAB")
#' format_taxon("Pseudomonas_B")
#' format_taxon("Pseudomonas")
format_taxon <- function(x) {
  x %>%
    stringr::str_trim() %>%
    { 
      dplyr::case_when(
        # Case 1: ALL uppercase and/or digits → leave unchanged
        stringr::str_detect(., "^[A-Z0-9]+$") ~ .,
        
        # Case 2: underscore + single letter suffix
        stringr::str_detect(., "_[A-Za-z]$") ~ 
          stringr::str_replace(., "^([^_]+)_([A-Za-z])$", "*\\1* \\2"),
        
        # Case 3: everything else → italicize whole thing
        TRUE ~ paste0("*", ., "*")
      )
    } %>%
    stringr::str_trim()
}
