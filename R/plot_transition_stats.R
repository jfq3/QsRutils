
#' Plot DADA2 Transition Stats
#' 
#' Extracts a QIIME2/DADA2 transition stats file and returns a plot
#'
#' @param trans_stats.qza The transitions stats file output by QIIME2 DADA2 
#' @details Beginning with QIIME2 version 2025.7 the DADA2 plugin requires outut of a compressed  (qza) file of the transition stats. This function makes a ggplot plot from the data in that file. It is useful in determining how well DADA2 has corrected read errors.
#' @returns A ggplot of the transition probabilities
#' @export
#'
#' @examples
#' \dontrun{}
plot_transition_stats <- function(trans_stats.qza) {
  
  # Find and extract Errorstats.tsv
  zip_list <- unzip(trans_stats.qza, list = TRUE)
  target <- zip_list$Name[grep("Errorstats\\.tsv$", zip_list$Name)]
  d <- unzip(trans_stats.qza, files = target)
  
  # Read file (skip QIIME2 metadata row)
  df <- readr::read_tsv(d, comment = "#", show_col_types = FALSE)
  
  # Clean + reshape forward reads
  df_f <- df |>
    dplyr::filter(F_from != F_to) |>  # keep only error transitions
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
  p <- ggplot2::ggplot(df_f, aes(x = Q)) +
    
    # Observed points
    geom_point(aes(y = observed),
               alpha = 0.5,
               size = 1.2) +
    
    # Estimated model (solid line)
    geom_line(aes(y = estimated),
              linewidth = 0.6) +
    
    # Nominal error (dashed line)
    geom_line(aes(y = nominal),
              linetype = "dashed",
              linewidth = 0.6,
              color = "red") +
    
    facet_wrap(~ transition, ncol = 3) +
    
    scale_y_log10(
      name = "Error rate (log scale)",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::label_scientific()
    ) +
    
    scale_x_continuous(
      name = "Quality score (Phred Q)"
    ) +
    
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      legend.position = "none",
    )
  
  p
}

