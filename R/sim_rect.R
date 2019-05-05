#' @title Simulation of postholes
#'
#' @description \code{sim_rect} allows the user to test the main function \code{find_rect()}
#' with simulated data. The number of houses as well as the number of rows and
#' crossrows can be choosen freely. The function allows furthermore for
#' uncertainty in measurements and placements of postholes and also simulates
#' the common situation that not all formerly existing postholes were
#' discovered. Finally, to add another layer of fuzziness, the function allows
#' to add further random points, because in real- life situations it is
#' virtually impossible to ascribe every posthole to a certain building.
#'
#' @param lower Number marking the lower left corner of the area of interest.
#' @param upper Number marking the lower left corner of the area of interest.
#' @param housenumber Number of simulated houses.
#' @param mindist Decimal, minimum distance between postholes.
#' @param maxdist Decimal, maximum distance between postholes.
#' @param pointjitter Decimal, defining the jittering of points (= postholes).
#' @param preserve Number, defining the percentage of points to be preserved
#'   (0-100).
#' @param maxrows Maximum number of rows and crossrows (minimum number is always
#'   2).
#' @param rp Number of additional random points (= postholes).
#' @param outputmethod Character string, defining the format of the output list.
#'   Options include shp (SpatialPointDataframe, SpatialLinesDataframe),
#'   spatstat, spatial
#'
#' @return A list of points and lines in various formats for display in R or
#'   QGIS.
#'
#' @examples
#' # Simulating 50 houses with 2-5 rows and crossrows where the posts are in a 
#' # distance of 2-5 units to each other in a window of 200 x 200 units. The
#' # points are jittered by 0.3 units, and an additional 1000 points are
#' # generated. Only 80% of all points are preserved. The output is a list
#' # spatstat-points and -lines. 
#' library(spatstat)
#' library(maptools)
#' library(sp)
#' library(data.table)
#' sim <- sim_rect(lower = 0, upper = 200, housenumber = 20,  maxrows = 5,
#' mindist = 2, maxdist = 5, pointjitter = 0.3, preserve = 80,
#' rp = 100, outputmethod = "spatstat") 
#' plot(sim$points)
#' plot(sim$lines, add=TRUE)
#'
#'
#' @importFrom stats runif
#' @importFrom grDevices chull
#' @importFrom maptools as.SpatialLines.psp
#' @importFrom maptools as.SpatialPoints.ppp
#' @importFrom sp SpatialLinesDataFrame
#' @importFrom sp SpatialPointsDataFrame
#' @importFrom spatstat as.psp
#' @importFrom spatstat owin
#' @importFrom spatstat ppp
#' @importFrom spatstat rjitter
#' @importFrom spatstat rthin
#'    
#' @export

sim_rect<-function(lower,upper,housenumber,maxrows,mindist,maxdist,pointjitter,preserve,rp,outputmethod) {
  houses <- NULL
  house_sim <- NULL
  rectangle_nr <- 0
  rec_nr <- NULL
  for (i in 1:housenumber) {
    house <- NULL
    rec_nr <- i
    # centroid of simulated house
    x<-sample((lower+maxrows*maxdist):(upper+maxrows*maxdist), size = 1)
    y<-sample((lower+maxrows*maxdist):(upper+maxrows*maxdist), size = 1)
    
    # actual distance between rows and crossrows per house
    rowdist <- sample(mindist:maxdist, size = 1)
    crossrowdist <- sample(mindist:maxdist, size = 1)
    
    # actual number of rows and crossrows per house
    if (maxrows < 3) {
      crossrows = 2
      rows = 2
    } else {
      crossrows<-sample(2:maxrows, size = 1)
      rows<-sample(2:maxrows, size = 1)
    }

    # cumulative sum of rectangles comprised by houses
    rectangle_nr<-(crossrows-1)*(rows-1)+rectangle_nr
        
    # defining the house orientation
    theta<-sample(1:180, size = 1)

    for (j in 1:crossrows) {
      for (m in 1:rows) {
        tempx<- -(-0.5*(rows+1)+m)*rowdist
        tempy<- -(-0.5*(crossrows+1)+j)*crossrowdist
        rotatedx <- tempx*cos(theta) - tempy*sin(theta);
        rotatedy <- tempx*sin(theta) + tempy*cos(theta);
        ax_rot<- rotatedx + x
        ay_rot<- rotatedy + y
        house<-rbind(house,data.frame(rec_nr,x=ax_rot,y=ay_rot))
      }
    }
    # defining the outline of the house
    houses<-rbind(houses,house)
    house_cps <- grDevices::chull(house$x,house$y)
    cps_nr <-length(house_cps)
    for (r in 1:cps_nr) {
      if (r < cps_nr) {
        r2=r+1
      } else {
        r2=1
      }
      house_sim <- rbind(house_sim,data.frame(marks=as.factor(rec_nr),x0=house[house_cps[r],]$x,
                                    y0=house[house_cps[r],]$y,x1=house[house_cps[r2],]$x,
                                    y1=house[house_cps[r2],]$y))
    }
  }
  
  # define upper limits for window for point and line patterns in spatstat
  upperlimits <- upper + 2 * maxrows * maxdist

  if (rp > 0) {
  # generate additional randompoints
  rand_x<-stats::runif(rp, lower, upperlimits) # x-coordinates
  rand_y<-stats::runif(rp, lower, upperlimits) # y-coordinates
  houses<-rbind(houses,data.frame(rec_nr=99999,x=rand_x,y=rand_y))
  }
  
  
  
  # create an object of type pointpattern in spatstat
  X <- spatstat::ppp(marks=as.factor(houses$rec_nr),houses$x, houses$y, 
                     window=spatstat::owin(c(lower,upperlimits),c(lower,upperlimits)))

  # create an object of type linepattern in spatstat
  X_line <- spatstat::as.psp(house_sim, window=spatstat::owin(c(lower,upperlimits),c(lower,upperlimits)))

if (preserve < 100) {
  # thin the dataset
  X_thinned <- spatstat::rthin(X,preserve/100)
  } else {
  X_thinned<-X
  }
  
  if (pointjitter > 0) {
  # jitter the points
  X_jittered <- spatstat::rjitter(X_thinned,pointjitter)
  } else {
  X_jittered<-X_thinned
  }
  
   #output <- cat(sep = '',"\n",c, " houses with ", rectangle_nr, " rectangles were simulated.\n\n")
  
  output<-NULL
 
  # Spatial
  sim_points <- maptools::as.SpatialPoints.ppp(X_jittered)  
  sim_lines <- maptools::as.SpatialLines.psp(X_line)
  
  # SpatialPointDataframe
  sim_points_shp <- sp::SpatialPointsDataFrame(sim_points,data.frame(id=1:length(sim_points)))
  
  # SpatialLinesDataframe
  sim_lines_shp <- sp::SpatialLinesDataFrame(sim_lines,data.frame(id=1:length(sim_lines)))


  # Define output method
  if (outputmethod=="spatstat") {
    points <- X_jittered
    lines <- X_line
  } else if (outputmethod=="Spatial") {
    points <- sim_points
    lines <- sim_lines
  } else if (outputmethod=="shape") {
      points <- sim_points_shp
     lines <- sim_lines_shp
   }
  
  return(list(points = points, lines = lines))
}
