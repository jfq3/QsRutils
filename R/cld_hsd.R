#' Make CLD tibble from Tukey HSD Results
#' 
#' Makes a tibble for adding compact letter assignments to a boxplot using HSD.test results.
#' 
#' @param hsd_rslt /The result of the HSD>test function of package agricolae
#' 
#' @return A tibble with columns for treatment groups (x), the 75th percintile of the treatment groups (y), and the CLD letters indicating significantly different treatments.
#' 
#' @details
#' hsd_rslt must be created with agricolae::HSD.test
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' data("iris")
#' model <- lm(Petal.Length ~ Species, data = iris)
#' hsd_rslt <- agricolae::HSD.test(model, trt="Species")
#' cld_hsd(hsd_rslt)
#' }
#' 
cld_hsd <- function(hsd_rslt) {
  tk_means <- hsd_rslt$means |> 
    tibble::rownames_to_column(var = "x") |> 
    dplyr::select(x, Q75) |> 
    dplyr::rename(y = "Q75")
  tk_grps <- hsd_rslt$groups |> 
    tibble::rownames_to_column(var = "x") |> 
    dplyr::select(x, groups)
  ltr_df <- dplyr::inner_join(tk_means, tk_grps, by = "x") 
  return(ltr_df)
}
