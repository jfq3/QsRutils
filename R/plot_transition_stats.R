
utils::globalVariables(c(
  "F_Estimated", "F_Nominal", "F_Observed", "F_Qual", "F_Transition",
  "F_from", "F_to",
  "R_Estimated", "R_Nominal", "R_Observed", "R_Qual", "R_Transition",
  "R_from", "R_to",
  "Q", "direction", "estimated", "nominal", "observed", "transition"
))

#' Plot DADA2 Transition Stats
#'
#' Extracts a QIIME2/DADA2 transition stats file and returns a plot.
#'
#' @param trans_stats.qza The transitions stats file output by QIIME2 DADA2
#' @details Beginning with QIIME2 version 2025.7 the DADA2 plugin requires
#'   output of a compressed (qza) file of the transition stats. This function
#'   makes a ggplot plot from the data in that file. It is useful in
#'   determining how well DADA2 has corrected read errors.
#' @returns A ggplot of the transition probabilities
#' @export
#' @importFrom utils unzip
#' @importFrom readr read_tsv
#' @importFrom scales trans_breaks label_scientific
#' @importFrom ggplot2 ggplot aes geom_point geom_line facet_wrap
#'   scale_y_log10 scale_x_continuous theme_bw theme element_text
#'   element_rect element_blank element_line
#'
#' @examples
#' \dontrun{
#'   plt <- plot_transition_stats("path/to/transition_stats.qza")
#'   print(plt)
#' }
plot_transition_stats <- function(trans_stats.qza) {

  # Find and extract Errorstats.tsv
  zip_list <- utils::unzip(trans_stats.qza, list = TRUE)
  target <- zip_list$Name[grep("Errorstats\\.tsv$", zip_list$Name)]
  d <- utils::unzip(trans_stats.qza, files = target)

  # Read file (skip QIIME2 metadata row)
  df <- readr::read_tsv(d, comment = "#", show_col_types = FALSE)

  # Clean + reshape forward reads
  df_f <- df |>
    dplyr::filter(F_from != F_to) |>
    dplyr::mutate(
      direction = "Forward",
      transition = F_Transition,
      Q = F_Qual,
      observed = F_Observed,
      estimated = F_Estimated,
      nominal = F_Nominal,
      .keep = "none"
    )

  # Clean + reshape reverse reads
  df_r <- df |>
    dplyr::filter(R_from != R_to) |>
    dplyr::mutate(
      direction = "Reverse",
      transition = R_Transition,
      Q = R_Qual,
      observed = R_Observed,
      estimated = R_Estimated,
      nominal = R_Nominal,
      .keep = "none"
    )

  # Combine
  plot_df <- dplyr::bind_rows(df_f, df_r) |>
    dplyr::mutate(
      transition = factor(transition),
      direction = factor(direction)
    )

  # ---- Plot ----
  ggplot2::ggplot(df_f, ggplot2::aes(x = Q)) +

    ggplot2::geom_point(ggplot2::aes(y = observed),
               alpha = 0.5,
               size = 1.2) +

    ggplot2::geom_line(ggplot2::aes(y = estimated),
              linewidth = 0.6) +

    ggplot2::geom_line(ggplot2::aes(y = nominal),
              linetype = "dashed",
              linewidth = 0.6,
              color = "red") +

    ggplot2::facet_wrap(~ transition, ncol = 3) +

    ggplot2::scale_y_log10(
      name = "Error rate (log scale)",
      breaks = scales::trans_breaks("log10", \(x) 10^x),
      labels = scales::label_scientific()
    ) +

    ggplot2::scale_x_continuous(
      name = "Quality score (Phred Q)"
    ) +

    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "white"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      legend.position = "none"
    )
}
