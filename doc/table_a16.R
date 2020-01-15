## ---- echo=FALSE--------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(QsRutils)))
suppressWarnings(suppressPackageStartupMessages(library(knitr)))

## -----------------------------------------------------------------------------
part1 <- matrix(data = NA, nrow = 10, ncol = 10)
for (i in 1:10) {
  for (j in 1:10) {
    x <- rad2deg(arc_sine(((i-1)/10)+((j-1)/100)))
    part1[i,j] <- formatC(x, digits = 2, format = "f")
  }
}
col1 <- seq(from = 0, to = 0.9, by = 0.1)
part1 <- data.frame(col1, part1)
colnames(part1) <- c("%", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")

part2 <- matrix(data = NA, nrow = 98, ncol = 10)
for (i in 1:98) {
  for (j in 1:10) {
    x <- rad2deg(arc_sine(i+((j-1)/10)))
    part2[i,j] <-formatC(x, digits = 2, format = "f")
  }
}
col1 <- seq(from = 1, to = 98, by = 1)
part2 <- data.frame(col1, part2)
colnames(part2) <- c("%", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")

part3 <- matrix(data = NA, nrow = 11, ncol = 10)
for (i in 1:11) {
  for (j in 1:10) {
    m <- 99 + ((i-1)/10)
    n <- ((j-1)/100)
    if((m+n) <= 100) {
      x <- rad2deg(arc_sine(m+n))
       part3[i,j] <- formatC(x, digits = 2, format = "f")
    }
  }
}
col1 <- round(seq(from = 99, to = 100, by = 0.1), 1)
part3 <- data.frame(col1, part3)
colnames(part3) <- c("%", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")

a16 <- rbind(part1, part2, part3)

knitr::kable(a16, align = "l", caption = "Table of Angles Corresponding to Percentages, after Bliss(1937).")

