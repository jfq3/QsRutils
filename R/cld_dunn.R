#' CLDs from DUNN Generate compact letter displays from Dunn.test results
#' 
#' @param dunn_rslt The result of the function dunn.test::dunn.test
#' @param significance The alpha level for  statistical significance
#' @return A character vector of CLDs - groups of significantly different treatment groups
#' @export
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom dplyr rowwise
#' @importFrom dplyr select
#' @importFrom stats as.dist
#' @importFrom stringr str_split_1
#' @importFrom stringr str_trim
#' @importFrom tibble tibble
#' @examples
#' \dontrun{
#'
#' }
#' 
cld_dunn <- function(dunn_rslt, significance = 0.05) {
  comparison <- comp1 <- comp2 <- p.adj <- NULL
  df_mat <- tibble(dunn_rslt$res$Comparison, dunn_rslt$res$P.adj)
  colnames(df_mat) <- c("comparison", "p.adj")
  df_mat <- df_mat %>% 
    rowwise() %>% 
    mutate(comp1=str_split_1(comparison, "-")[1],
           comp1=str_trim(comp1),
           comp2=str_split_1(comparison, "-")[2],
           comp2=str_trim(comp2)) %>% 
    select(comp1, comp2, p.adj) %>% 
    arrange(comp1)
  trts <- unique(c(df_mat$comp1, df_mat$comp2))
  stat <- df_mat$p.adj
  stat.d <- matrix(NA, nrow = length(trts), ncol = length(trts))
  k=1
  for (j in 1:(nrow(stat.d)-1)) {
    for (i in ((j+1):nrow(stat.d))) {
      stat.d[i,j] <- stat[k]
      k=k+1
    }
  }
  stat.d <- as.matrix(stat.d)
  rownames(stat.d) <- colnames(stat.d) <- trts
  stat.d <- as.dist(stat.d)
  letter.assignments <- multcompView::multcompLetters(stat.d, threshold = significance)
  return(list(p_adj_matrix = stat.d, clds = letter.assignments$Letters))
}
