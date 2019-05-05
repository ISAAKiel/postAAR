#' @title Expected number of rectangles
#'
#' @description \code{exp_rect} allows the user to determine the expected
#' number of rectangles according to the formula by Fletcher and Lock 1984. This
#' is important because time of computation rises exponentially with the number
#' of rectangles. With more than 4000-5000 expected rectangles, computation 
#' will take unreasonable long (or more likely trigger a memory error).
#' 
#' In comparison to the original formula which had only "0.5" in the bottom of
#' the formula leading to "NRE", we found that this tends to doubling the true
#' expected number. Therefore, we changed that this term to "0.25" which produces
#' numbers in accordance with simulated runs.
#' 
#' @param x_coord Vector of x-coordinates.
#' @param y_coord Vector of y-coordinates.
#' @param mindist Decimal, minimum distance between postholes.
#' @param maxdist Decimal, maximum distance between postholes.
#' @param pointjitter Decimal, defining the jittering of points.
#'
#' @return A string explaining the number of expected rectangles.
#'
#' @examples
#' library(spatstat)
#' library(maptools)
#' library(sp)
#' library(data.table)
#' sim <- sim_rect(lower = 0, upper = 200, housenumber = 0,  maxrows = 5,
#' mindist = 2, maxdist = 5, pointjitter = 0.3, preserve = 100,
#' rp = 16200, outputmethod = "spatstat")
#' 
#' # exporting the coordinates
#' sim_coords <- spatstat::coords(sim$points)
#' x_coord <- sim_coords$x
#' y_coord <- sim_coords$y
#' 
#' # calculating the number of expected rectangles
#' exp_rect(x_coord, y_coord, mindist = 2, maxdist = 5, pointjitter = 0.3)
#' 
#' @export
  exp_rect <- function(x_coord, y_coord, mindist, maxdist, pointjitter) {
  points <- cbind(x_coord, y_coord) 
  hull <- geometry::convhulln(points, options = 'FA') # convex hull of points
  A_hull <- hull$vol # area of convex hull
  pn <- length(x_coord) # number of points/postholes
  
  # formula according to Fletcher and Lock
  RM <- (maxdist * maxdist - mindist * mindist) * pi * pn * (pn -1) / A_hull
  RMU <- (pn - 2) * 2 * (maxdist - mindist) * pointjitter / A_hull
  RLA <- -2.41 * pointjitter * pointjitter * (pn - 2) / A_hull
  SI <- RMU * (1 - exp(RLA))
  NRE <- round(0.25 * RM * exp(-SI) * SI / (1 - SI)) # NOTE: the original term "0.5" is replaced here with "0.25"
  return(cat(sep = '',"\n","With the given parameters (minimal distance: ", 
             mindist, ", maximal distance: ", maxdist, ", tolerance: ", pointjitter, 
              "), the expected number of rectangles is:\n\n", NRE,".\n"))
}
  