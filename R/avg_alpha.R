#' Average Alpha Diversity (faster implementation)
#'
#' Calculates alpha-diversity metrics from n samplings of an OTU table to a constant number of counts per sample.
#'
#' This implementation focuses on speed:  
#' - vectorizes the alpha-metric calculations (avoids repeated calls to vegan::diversity and vegan::specnumber)  
#' - uses base R operations (rowSums, logical ops, arithmetic) which are much faster in tight loops  
#' - optionally parallelizes replicates across cores (platform-aware).  
#'
#' @param otu An OTU table as a data frame or matrix with samples as rows and taxa as columns.
#' @param sampling_depth The number of counts per sample in the sampled OTU table
#' @param iterations The number of times the OTU table should be sampled.
#' @param sum_method Method ("median" or "mean") for summarizing replication results.
#' @param ncores Number of cores to use for parallel execution. Default 1 (no parallelism).
#'
#' @details
#' The OTU data frame supplied must be in typical vegan format: samples as row names and taxa as column names.
#' The minimum row sum must be greater than or equal to the sampling depth.
#' 
#' By default the sum_method is mean. For a similar function in QIIME2, the default sum_method is median.
#'
#' @return Returns a dataframe with Shannon, Observed, Pielou, Simpson and Inverse Simpson for each sample in an OTU table.
#' @export
#'
#' @examples
#' {
#' data(BCI, package = "vegan")
#' otu <- BCI[rowSums(BCI) > 400, ]
#' avg_alpha(otu, sampling_depth = 400, iterations = 100)
#' }
avg_alpha <- function(otu, sampling_depth, iterations = 100, sum_method = c("median", "mean"), ncores = 1) {
  sum_method <- match.arg(sum_method)
  summary_fun <- if (sum_method == "median") stats::median else base::mean

  if (!is.matrix(otu) && !is.data.frame(otu)) stop("otu must be a matrix or data.frame")
  otu <- as.matrix(otu)
  if (any(is.na(otu))) stop("otu contains NA values; please remove or impute them first")

  if (nrow(otu) == 0L) return(data.frame())

  if (min(rowSums(otu)) < sampling_depth) {
    stop("All samples must have at least 'sampling_depth' counts (min row sum < sampling_depth).")
  }

  a <- nrow(otu)
  sample_names <- rownames(otu)
  metric_names <- c("Shannon", "Observed", "Pielou", "Simpson", "InvSimpson")

  # Single-iteration worker (vectorized per-iteration).
  # After rrarefy(), every row sum equals sampling_depth, so we divide by the
  # scalar rather than recomputing rowSums() — one fewer full matrix pass.
  compute_once <- function(i) {
    otu_r <- rrarefy_cpp(otu, sampling_depth)

    p <- otu_r / sampling_depth  # faster than otu_r / rowSums(otu_r)

    # Shannon: -sum(p * log(p)), 0*log(0) defined as 0.
    # Avoid a second matrix allocation by reusing p: temporarily set zeros to 1
    # so log(1) = 0, then the product p * log(p_safe) is still 0 at those cells.
    p_safe <- p
    nz <- p > 0
    p_safe[!nz] <- 1
    shannon <- -rowSums(p * log(p_safe))

    observed  <- rowSums(otu_r > 0L)

    # Pielou's evenness: Shannon / log(S); undefined when S <= 1
    pielou <- rep(NA_real_, a)
    idx <- observed > 1L
    pielou[idx] <- shannon[idx] / log(observed[idx])

    simpson    <- rowSums(p * p)
    invsimpson <- 1 / simpson

    cbind(shannon, observed, pielou, simpson, invsimpson)
  }

  if (!is.numeric(ncores) || length(ncores) != 1 || ncores < 1) ncores <- 1L
  ncores <- as.integer(ncores)

  if (ncores == 1L) {
    if (sum_method == "mean") {
      # Accumulate directly — avoids allocating the full samples x metrics x iterations array.
      acc <- matrix(0, nrow = a, ncol = length(metric_names))
      for (i in seq_len(iterations)) acc <- acc + compute_once(i)
      result_mat <- acc / iterations
    } else {
      # Median requires all replicates in memory.
      reps <- lapply(seq_len(iterations), compute_once)
      arr  <- array(unlist(reps), dim = c(a, length(metric_names), iterations))
      result_mat <- apply(arr, c(1, 2), stats::median)
    }
  } else {
    # Parallel path — workers need vegan, the otu matrix, sampling_depth, and
    # compute_once exported explicitly (PSOCK workers don't inherit the closure).
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makePSOCKcluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, c("otu", "sampling_depth", "a", "compute_once"),
                              envir = environment())
      reps <- parallel::parLapply(cl, seq_len(iterations), compute_once)
    } else {
      reps <- parallel::mclapply(seq_len(iterations), compute_once, mc.cores = ncores)
    }

    arr        <- array(unlist(reps), dim = c(a, length(metric_names), iterations))
    result_mat <- apply(arr, c(1, 2), summary_fun)
  }

  result_df <- as.data.frame(result_mat, stringsAsFactors = FALSE)
  colnames(result_df) <- metric_names
  rownames(result_df) <- sample_names

  return(result_df)
}
