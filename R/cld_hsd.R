#' Make CLD tibble from Tukey HSD Results
#' 
#' Makes a tibble for adding compact letter assignments to a boxplot using HSD.test results.
#' 
#' @param hsd_rslt The result of the HSD>test function of package agricolae
#' @param y_pos The y-position in relation to the boxplots. Choices are at the top of the box ("boxtop", the default) or at the maximum group value ("max").
#' @return A tibble with columns for treatment groups (x), the y-positions of the treatment CLD (y), and the CLD letters indicating significantly different treatments.
#' 
#' @details
#' hsd_rslt must be created with agricolae::HSD.test
#' 
#' @export
#' 
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom dplyr inner_join
#' 
#' @examples
#' \dontrun{
#' data("iris")
#' model <- lm(Petal.Length ~ Species, data = iris)
#' hsd_rslt <- agricolae::HSD.test(model, trt="Species")
#' cld_hsd(hsd_rslt)
#' }
#' 
cld_hsd <- function(hsd_rslt, y_pos = "boxtop") {
  x <- Q75 <- groups <- NULL
  if (y_pos == "boxtop") {
    tk_means <- hsd_rslt$means |> 
    tibble::rownames_to_column(var = "x") |> 
    dplyr::select(x, Q75) |> 
    dplyr::rename(y = "Q75")
  } else {
    tk_means <- hsd_rslt$means |> 
      tibble::rownames_to_column(var = "x") |> 
      dplyr::select(x, Max) |> 
      dplyr::rename(y = "Max")
  }

  tk_grps <- hsd_rslt$groups |> 
    tibble::rownames_to_column(var = "x") |> 
    dplyr::select(x, groups)
  ltr_df <- dplyr::inner_join(tk_means, tk_grps, by = "x") 
  return(ltr_df)
}
