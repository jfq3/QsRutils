#' Averaged Sub-sampled Dissimilarity Matrices
#' 
#' The function computes the dissimilarity matrix of a dataset
#'   multiple times using \code{\link{vegdist}} while randomly subsampling the
#'   dataset each time. All of the subsampled iterations are then averaged
#'   (mean) to provide a distance matrix that represents the average of multiple
#'   subsampling iterations. This emulates the behavior of the distance matrix
#'   calculator within the Mothur microbial ecology toolkit.
#' @usage avg_dist(x, sample, distfun = vegdist, meanfun = mean,
#' transf = NULL, iterations = 100, dmethod = "bray",
#' diag = TRUE, upper = TRUE, ncores = 1L, ...)
#' @param x Community data matrix
#' @param sample The subsampling depth to be used in each iteration. Samples that do not meet this threshold will be removed from the analysis, and their identity returned to the user in stdout.
#' @param distfun The dissmilarity matrix function to he used. Default is the vegan vegdist
#' @param meanfun The calculation to use for the average (mean or median)
#' @param transf Transformation function to apply to the distance matrix. Default is NULL.
#' @param iterations Number of subsampling iterations. Default is 100.
#' @param dmethod Distance method. Default is "bray".
#' @param diag Whether to include the diagonal of the distance matrix. Default is TRUE.
#' @param upper Whether to include the upper triangle of the distance matrix. Default is TRUE.
#' @param ncores Number of cores to use for parallel computation via
#'   \code{\link[parallel]{makeCluster}}. Default is \code{1L} (no parallelism).
#'   Values greater than 1 spawn a socket cluster, which works on all platforms
#'   including Windows.
#' @param ... Any additional arguments to add to the distance function of mean/median function specified.
#' @return An averaged distance matrix.
#' @details
#' **Parallelism.** Each iteration is independent, so the work can be
#' distributed across cores with \code{ncores > 1}. Whether this is
#' worthwhile depends on the cost of each iteration relative to the fixed
#' overhead of spawning a socket cluster (typically 1–2 seconds).
#'
#' Benchmarks on a Windows machine with 100 samples × 5000 OTUs
#' (Bray-Curtis, \code{sample = 20000}) gave the following observed speedups
#' relative to \code{ncores = 1}:
#'
#' \tabular{rrr}{
#'   \strong{ncores} \tab \strong{50 iters} \tab \strong{200 iters} \cr
#'   1  \tab 1.00× (48 s)  \tab 1.00× (204 s) \cr
#'   2  \tab 1.44× (33 s)  \tab 1.48× (138 s) \cr
#'   4  \tab 1.91× (25 s)  \tab 2.04× (100 s) \cr
#'   8  \tab 2.12× (23 s)  \tab 2.37× (86 s)  \cr
#' }
#'
#' Observed speedups are well below the theoretical maximum due to
#' inter-process communication overhead in socket clusters (the only cluster
#' type available on Windows). The relative benefit of parallelism improves
#' with more iterations because the fixed cluster overhead is amortised across
#' a larger amount of work. For small problems (e.g., <50 samples and
#' \code{iterations} ≤ 100 with the default \code{vegdist}), the per-iteration
#' cost is low enough that \code{ncores = 1L} is likely faster.
#' @export
#' @author Geoffrey Hannigan, with some minor tweaks by Gavin L. Simpson.
#' @note The function builds on the function \code{\link{rrarefy}} and an additional distance matrix function (e.g. \code{\link{vegdist}}) to add more meaningful representations of distances among randomly subsampled datasets by presenting the average of multiple random iterations. This functionality has been utilized in the Mothur standalone microbial ecology toolkit, see https://mothur.org/wiki/Dist.shared.
#' @seealso [vegdist()]
#' @seealso [rarrefy()]
#' @keywords multivariate
#' @examples
#' # Import an example count dataset
#' data(BCI)
#' # Test the base functionality
#' mean.avg.dist <- avg_dist(BCI, sample = 50, iterations = 10)
#' # Test the transformation function
#' mean.avg.dist.t <- avg_dist(BCI, sample = 50, iterations = 10, transf = sqrt)
#' # Test the median functionality
#' median.avg.dist <- avg_dist(BCI, sample = 50, iterations = 10, meanfun = median)
#' # Print the resulting tables
#' head(as.matrix(mean.avg.dist))
#' head(as.matrix(mean.avg.dist.t))
#' head(as.matrix(median.avg.dist))
#' # Run example to illustrate low variance of mean, median, and stdev results
#' # Mean and median std dev are around 0.05
#' sdd <- avg_dist(BCI, sample = 50, iterations = 100, meanfun = sd)
#' summary(mean.avg.dist)
#' summary(median.avg.dist)
#' summary(sdd)
#' # Test for when subsampling depth excludes some samples
#' # Return samples that are removed for not meeting depth filter
#' depth.avg.dist <- avg_dist(BCI, sample = 450, iterations = 10)
#' # Print the result
#' depth.avg.dist

avg_dist <- function (x, sample, distfun = vegdist, meanfun = mean, transf = NULL, 
                     iterations = 100, dmethod = "bray", diag = TRUE, upper = TRUE,
                     ncores = 1L, ...) 
{
  if (missing(sample)) {
    stop("Subsampling depth must be supplied via argument 'sample'")
  }
  else {
    if (!(is.numeric(sample) && sample > 0L)) {
      stop("Invalid subsampling depth; 'sample' must be positive & numeric")
    }
  }
  if (!is.numeric(iterations)) {
    stop("Invalid iteration count; must be numeric")
  }
  inputcast <- x
  distfun <- match.fun(distfun)
  if (!is.null(transf)) {
    transf <- match.fun(transf)
  }
  minobs <- min(x[x > 0])
  if (minobs > 1) 
    warning(gettextf("most observed count data have counts 1, but smallest count is %d", 
                     minobs))

  iter_fun <- function(i) {
    mat <- suppressWarnings(rrarefy(inputcast, sample = sample))
    mat <- mat[rowSums(mat) >= sample, , drop = FALSE]
    if (!is.null(transf)) {
      mat <- transf(mat)
    }
    as.matrix(distfun(mat, method = dmethod, diag = TRUE, upper = TRUE, ...))
  }

  ncores <- as.integer(ncores)
  if (ncores > 1L) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, varlist = c("inputcast", "sample", "transf",
                                            "distfun", "dmethod"),
                            envir = environment())
    parallel::clusterEvalQ(cl, library(vegan))
    distlist <- parallel::parLapply(cl, seq_len(iterations), iter_fun)
  } else {
    distlist <- lapply(seq_len(iterations), iter_fun)
  }

  meanfun <- match.fun(meanfun)
  rnames  <- row.names(distlist[[1]])
  d       <- dim(distlist[[1]])
  afunc   <- array(unlist(distlist), c(d, length(distlist)))

  # rowMeans / rowMedians are far faster than apply() for the common cases
  if (identical(meanfun, mean)) {
    output <- rowMeans(afunc, dims = 2L)
  } else {
    output <- apply(afunc, 1:2, meanfun, ...)
  }

  colnames(output) <- rownames(output) <- rnames
  dropsamples <- setdiff(row.names(x), rnames)
  if (length(dropsamples) > 0L) {
    warning(gettextf("The following sampling units were removed because they were below sampling depth: %s", 
                     paste(dropsamples, collapse = ", ")))
  }
  output <- as.dist(output, diag = diag, upper = upper)
  attr(output, "call") <- match.call()
  attr(output, "method") <- paste("avg_dist", dmethod)
  return(output)
}