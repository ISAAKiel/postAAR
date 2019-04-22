##Archaeological analysis=group
##Line_layer=vector
##MBR_hull=output vector

library("spatstat")
library("maptools")
library("sp")
library("data.table")
library("geometry")

psp=as(as(Line_layer, "SpatialLines"),"psp")
psp_coords<-coords(psp)

MBR2 <- function(points) {
  tryCatch({
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
    return(as.data.frame(cbind(x[c(1,2,2,1,1),k], y[c(1,1,2,2,1),k]) %*% rbind(v[k,], w[k,])))
  }, error = function(e) {
    assign('points', points, .GlobalEnv)
    stop(e)  
  })
}

hull<-MBR2(ppp_coords)
MBR_hull<-SpatialPointsDataFrame(hull,data.frame(id=1:nrow(hull)))