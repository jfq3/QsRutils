#' log_arc_sine
#'
#' Log of the arc-sine Transfromation of a Percentage
#'
#' @param x A percentage.
#'
#' @return The common logarithm of the arcsine transformation of x.
#' @export
#'
#' @examples
#'
#' log_arc_sine(x = 30.1)
#'
log_arc_sine <- function(x) {
  y <- log((asin(sqrt(x/100)))+1)
  return(y)
}

#' sqrt_arc_sine
#'
#' Square Root of the arc-sine of a Percentage
#'
#' @param x A percentage.
#'
#' @return The square root of the arcsine transformation of x.
#' @export
#'
#' @examples
#'
#' sqrt_arc_sine(30.1)
#'
sqrt_arc_sine <- function(x) {
  y <- sqrt(asin(sqrt(x/100)))
  return(y)
}

#' arc_sine
#'
#' Arcsine of a Percentage
#'
#' @param x A percentage.
#'
#' @return The arcsine transformation of x.
#' @export
#'
#' @examples
#'
#' arc_sine(30.1)
#'
arc_sine <- function(x) {
  y <- asin(sqrt(x/100))
  return(y)
}

#' rad2deg
#'
#' Radians to degrees
#'
#' @param x Angle in radians
#'
#' @return Angle in degrees.
#' @export
#' @examples
#'
#' rad2deg(pi * 0.5)
#'
rad2deg <- function(x) {
  y <- (x*180/pi)
  return(y)
}

#' deg2rad
#'
#' Degrees to radians
#'
#' @param x Angle in degrees
#' @return Angle in radians.
#' @export
#' @examples
#'
#' deg2rad(90)
#'
deg2rad <- function(x) {
  y <- x*(pi/180)
  return(y)
}

#'
