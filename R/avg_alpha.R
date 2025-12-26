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
#' data(BCI)
#' otu <- BCI[rowSums(BCI) > 400, ]
#' avg_alpha(otu, sampling_depth = 400, iterations = 100)
#' }
avg_alpha <- function(otu, sampling_depth, iterations = 100, sum_method = c("median", "mean"), ncores = 1) {
  sum_method <- match.arg(sum_method)
  summary_fun <- if (sum_method == "median") stats::median else base::mean
  
  if (!is.matrix(otu) && !is.data.frame(otu)) stop("otu must be a matrix or data.frame")
  otu <- as.matrix(otu)
  if (any(is.na(otu))) stop("otu contains NA values; please remove or impute them first")
  
  if (nrow(otu) == 0L) return(data.frame())  # nothing to do
  
  if (min(rowSums(otu)) < sampling_depth) {
    stop("All samples must have at least 'sampling_depth' counts (min row sum < sampling_depth).")
  }
  
  a <- nrow(otu)
  sample_names <- rownames(otu)
  metric_names <- c("Shannon", "Observed", "Pielou", "Simpson", "InvSimpson")
  
  # single-iteration worker (vectorized per-iteration)
  compute_once <- function(i) {
    otu_r <- vegan::rrarefy(otu, sampling_depth)
    # row sums should be equal to sampling_depth, but compute to be safe
    rs <- rowSums(otu_r)
    # probabilities
    p <- otu_r / rs
    
    # Shannon: -sum(p * log(p)) with p==0 contributing 0
    nz <- p > 0
    p_log_p <- matrix(0, nrow = nrow(p), ncol = ncol(p))
    p_log_p[nz] <- p[nz] * log(p[nz])
    shannon <- -rowSums(p_log_p)
    
    # Observed (richness) -- faster than vegan::specnumber
    observed <- rowSums(otu_r > 0)
    
    # Pielou's evenness: Shannon / log(S)
    pielou <- rep(NA_real_, length(observed))
    idx <- observed > 1
    pielou[idx] <- shannon[idx] / log(observed[idx])
    
    # Simpson: sum(p^2)
    p2 <- p * p
    simpson <- rowSums(p2)
    
    # InvSimpson: 1 / sum(p^2)
    invsimpson <- 1 / simpson
    
    cbind(shannon, observed, pielou, simpson, invsimpson)
  }
  
  # run iterations, optionally in parallel
  if (!is.numeric(ncores) || length(ncores) != 1 || ncores < 1) ncores <- 1L
  ncores <- as.integer(ncores)
  
  if (ncores == 1L) {
    reps <- lapply(seq_len(iterations), compute_once)
  } else {
    # Use parallel::mclapply on non-Windows; on Windows fall back to PSOCK cluster
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makePSOCKcluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      reps <- parallel::parLapply(cl, seq_len(iterations), function(i) {
        # import vegan in worker
        evalq({
          # library(vegan)
        }, envir = .GlobalEnv)
        compute_once(i)
      })
    } else {
      reps <- parallel::mclapply(seq_len(iterations), compute_once, mc.cores = ncores)
    }
  }
  
  # reps is a list of matrices (samples x metrics). Stack into an array: samples x metrics x iterations
  arr <- array(unlist(reps), dim = c(a, length(metric_names), iterations))
  dimnames(arr) <- list(sample_names, metric_names, paste0("rep", seq_len(iterations)))
  
  # summarize across iterations by sample and metric
  result_mat <- apply(arr, c(1, 2), summary_fun)
  
  result_df <- as.data.frame(result_mat, stringsAsFactors = FALSE)
  colnames(result_df) <- metric_names
  rownames(result_df) <- sample_names
  
  return(result_df)
}
