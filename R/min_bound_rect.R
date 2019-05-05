#' @title Minimum bounding rectangle
#' 
#' @description this function finds bounding rectangle. The main body of the 
#' function derives from this discussion on stackoverflow (dd). Credits go to blutt
#'
#' @param ppp_points   Points of spatstat::class ppp with marks as factor-level.
#' @param outputmethod Character string, defining the format of the output list.
#'   Options include shp (SpatialPointDataframe, SpatialLinesDataframe),
#'   spatstat, spatial
#'
#' @return
#' A list of lines delineating the minimum bounding rectangles of points
#' marked as belonging together. The format depends on the output method
#' chosen by the user.
#' 
#' 
#'
#' @examples
#' library(spatstat)
#' library(maptools)
#' library(sp)
#' library(data.table)
#' 
#' # setting up the example data set
#' sim <- sim_rect(lower = 0, upper = 200, housenumber = 20,  maxrows = 5,
#' mindist = 2, maxdist = 5, pointjitter = 0.3, preserve = 80,
#' rp = 200, outputmethod = "spatstat")
#' 
#' # exporting the coordinates
#' psp_coords <- spatstat::coords(sim$points)
#' 
#' # finding the rectangles
#' psp <- find_rect(psp_coords$x,psp_coords$y, 2, 6, 0.3, outputmethod = "shape")
#' 
#' # preparing the points
#' ppp_points <- as.ppp(as(psp$lines, "SpatialPointsDataFrame"))
#' 
#' # adding the marks
#' ppp_points$marks <- as.factor(ppp_points$marks$group_nr)
#' 
#' # performing the analysis
#' mbr <- min_bound_rect(ppp_points,outputmethod="spatstat")
#' plot(mbr)
#' 
#' @importFrom geometry convhulln
#' @importFrom sp Line
#' @importFrom sp Lines
#' @importFrom sp SpatialPointsDataFrame
#' @importFrom sp SpatialLinesDataFrame
#' @importFrom sp SpatialLines
#' @importFrom spatstat as.ppp
#' @importFrom spatstat as.psp
#' @importFrom spatstat coords
#' 
#' @export
#' 
min_bound_rect <- function(ppp_points, outputmethod) {
    marks_levels <- levels(ppp_points$marks)
    hulls <- NULL
    marks <- NULL
    for (t in marks_levels) {
      points_marks <- subset(ppp_points,marks==t)
      points <- spatstat::coords(points_marks)
      a2 <- geometry::convhulln(points, options = 'FA')
      
      e <- points[a2$hull[,2],] - points[a2$hull[,1],]            # Edge directions
      norms <- apply(e, 1, function(x) sqrt(x %*% x)) # Edge lengths
      
      v <- diag(1/norms) %*% as.matrix(e)                        # Unit edge directions
      
      
      w <- cbind(-v[,2], v[,1])                       # Normal directions to the edges
      
      # Find the MBR
      vertices <- as.matrix((points) [a2$hull, 1:2])    # Convex hull vertices
      minmax <- function(x) c(min(x), max(x))         # Computes min and max
      x <- apply(vertices %*% t(v), 2, minmax)        # Extremes along edges
      y <- apply(vertices %*% t(w), 2, minmax)        # Extremes normal to edges
      areas <- (y[1,]-y[2,])*(x[1,]-x[2,])            # Areas
      k <- which.min(areas)                           # Index of the best edge (smallest area)
      
      # Form a rectangle from the extremes of the best edge
      hull_line<-sp::Lines(list(sp::Line(cbind(x[c(1,2,2,1,1),k], y[c(1,1,2,2,1),k]) %*% rbind(v[k,], w[k,]))), ID = t)
      if (t==1) {
      hulls <- (list(hull_line))
      } else {
        hulls <- append(hulls, hull_line)
      }
    }
    hulls <- sp::SpatialLines(hulls)
      if (outputmethod == "Spatial") {
      return(hulls)
    } else if (outputmethod == "shape") {
      return(sp::SpatialLinesDataFrame(hulls, data.frame(id=1:length(hulls))), match.ID=FALSE)
    } else if (outputmethod == "spatstat") {
      return(spatstat::as.psp(hulls))
    }
}
