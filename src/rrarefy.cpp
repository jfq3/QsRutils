#include <Rcpp.h>
using namespace Rcpp;

//' Fast rarefaction via sequential hypergeometric sampling
//'
//' For each row, taxa counts are subsampled without replacement to exactly
//' \code{depth} total counts using the sequential hypergeometric algorithm:
//' iterate over taxa and draw from Hypergeometric(m = count_j,
//' n = remaining_pool - count_j, k = remaining_depth), updating both
//' running totals after each taxon. This is O(rows * cols) and produces
//' the same distribution as drawing \code{depth} items uniformly at random
//' without replacement from the pool of individuals.
//'
//' @param otu Numeric matrix of OTU counts (samples x taxa). Row sums must
//'   be >= \code{depth}.
//' @param depth Integer rarefaction depth.
//' @return Integer matrix of rarefied counts with the same dimensions and
//'   dimnames as \code{otu}.
//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix rrarefy_cpp(const NumericMatrix& otu, int depth) {
  int nr = otu.nrow();
  int nc = otu.ncol();
  IntegerMatrix out(nr, nc);

  for (int i = 0; i < nr; i++) {
    // Compute integer row sum (the total pool size for this sample)
    int pool = 0;
    for (int j = 0; j < nc; j++) pool += (int)otu(i, j);

    int rem_depth = depth; // draws still to assign

    for (int j = 0; j < nc; j++) {
      int c = (int)otu(i, j);

      if (c == 0 || rem_depth == 0) {
        out(i, j) = 0;
      } else if (pool == c) {
        // Only this taxon remains in the pool; assign all leftover draws to it.
        out(i, j) = rem_depth;
        rem_depth = 0;
      } else if (pool == rem_depth) {
        // Must take every remaining individual (no randomness left).
        out(i, j) = c;
        rem_depth -= c;
      } else {
        // Draw from Hypergeometric(white = c, black = pool - c, k = rem_depth).
        // R::rhyper(nn, mm, kk): nn = white balls, mm = black balls, kk = draws.
        int x = (int)R::rhyper((double)c, (double)(pool - c), (double)rem_depth);
        out(i, j) = x;
        rem_depth -= x;
      }

      pool -= c; // remove this taxon from the conceptual pool
    }
  }

  // Preserve row and column names
  out.attr("dimnames") = otu.attr("dimnames");
  return out;
}
