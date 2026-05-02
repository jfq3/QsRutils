#' An Experiment Level phyloseq Object
#'
#' Based on ITS2 sequences amplified from corn roots.
#'
#' @format A phyloseq object with otu_table, sample_data and tax_table. The sample_data variables are:
#' \describe{
#'   \item{P}{Phosporous level, H or L}
#'   \item{Genotype}{ One of three: 2, 3, and C}
#'   \item{Label}{A code for treatments: 2HR, 2LR, 3HR, 3LR, CHR, CLR}
#' }
"its.root"

#' A Data File in Long Format
#'
#' Used in Case 3 of the vignette make_comparisons
#'
#' @format A data file in long format used for a ggplot. The sample_data variables are:
#' \describe{
#'   \item{Treatment}{A code for genotype (2, 3, or C), P level (H or L) and sample type (R)}
#'   \item{Family}{One of the families in Gigasporaceae}
#'   \item{Percent}{Percent of total counts for family and treatment combination.}
#' }
"plot_df"


