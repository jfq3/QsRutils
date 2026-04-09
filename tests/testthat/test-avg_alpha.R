library(vegan)

test_that("rrarefy_cpp produces correct row sums", {
  data(BCI)
  otu <- as.matrix(BCI[1:5, ])
  result <- rrarefy_cpp(otu, 200)
  expect_equal(unname(rowSums(result)), rep(200L, 5))
})

test_that("rrarefy_cpp preserves dimnames", {
  data(BCI)
  otu <- as.matrix(BCI[1:3, 1:10])
  result <- rrarefy_cpp(otu, 200)
  expect_equal(dimnames(result), dimnames(otu))
})

test_that("rrarefy_cpp returns non-negative integers", {
  data(BCI)
  otu <- as.matrix(BCI[1:5, ])
  result <- rrarefy_cpp(otu, 200)
  expect_true(is.integer(result))
  expect_true(all(result >= 0L))
})

test_that("rrarefy_cpp mean richness matches vegan::rrarefy closely", {
  set.seed(6712)
  data(BCI)
  otu <- as.matrix(BCI[1:10, ])
  depth <- 200

  r_cpp   <- mean(replicate(300, rowSums(rrarefy_cpp(otu, depth) > 0)))
  r_vegan <- mean(replicate(300, rowSums(vegan::rrarefy(otu, depth) > 0)))

  # Means should agree to within ~1 species
  expect_lt(abs(r_cpp - r_vegan), 1.0)
})

test_that("avg_alpha returns correct dimensions and column names", {
  data(BCI)
  otu <- BCI[rowSums(BCI) > 400, ]
  result <- avg_alpha(otu, sampling_depth = 400, iterations = 10)

  expect_equal(nrow(result), nrow(otu))
  expect_equal(colnames(result), c("Shannon", "Observed", "Pielou", "Simpson", "InvSimpson"))
})

test_that("avg_alpha errors on insufficient row sums", {
  data(BCI)
  expect_snapshot(error = TRUE, avg_alpha(BCI, sampling_depth = 1e9, iterations = 5))
})
