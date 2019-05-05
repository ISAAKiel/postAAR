#' @title Find Rectangles
#'
#' @param a Vector of x-coordinates.
#' @param b Vector of y-coordinates.
#' @param x Decimal, minimum distance between postholes.
#' @param y Decimal, maximum distance between postholes.
#' @param z Decimal, defining the tolerance of points forming rectangles.
#' @param outputmethod Character string, defining the format of the output list.
#'   Options include shp (SpatialPointDataframe, SpatialLinesDataframe),
#'   spatstat, spatial
#'
#' @return A list of line and point spatstat objects.
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
#' # performing the search for rectangles
#' psp <- find_rect(psp_coords$x,psp_coords$y,2,5,0.3, outputmethod = "spatstat")
#' 
#' # plot the result
#' plot(psp$points)
#' plot(psp$lines, add=TRUE)

#' @useDynLib postAAR, .registration = TRUE
#' @importFrom Rcpp sourceCpp

#' @importFrom data.table rbindlist
#' @importFrom grDevices chull
#' @importFrom spatstat as.psp
#' @importFrom spatstat owin
#' @importFrom spatstat ppp
#' 
#' @export
#' 
find_rect <- function(a, b, x, y, z, outputmethod) {
  if (!is.loaded('postaar')) {
    dyn.load("./src/postAAR.so")
  }
  NH = 0
  NR = 0
  amin <- min(a)
  bmin <- min(b)
  a_min <- a - amin
  b_min <- b - bmin
  rec_nr <- NULL
  rectangles <- NULL
  
  a_ <- as.integer(round(a_min*100))
  b_ <- as.integer(round(b_min*100))
  z_ <- as.integer(round(z*100))
  x_ <- as.integer(round(x*100))
  y_ <- as.integer(round(y*100))
  NHO = length(a_)
  NSO = 1000000
  retvals<-.Fortran("postaar",
                    SMIN=as.integer(x_),
                    SMAX=as.integer(y_),
                    TOL=as.integer(z_),
                    NH=as.integer(NH),
                    NR=as.integer(NR),
                    NHO=as.integer(NHO),
                    NSO=as.integer(NSO),
                    IX=as.integer(a_),
                    IY=as.integer(b_),
                    POSTS=integer(NHO),
                    PX1=integer(NHO),
                    PY1=integer(NHO),
                    PX2=integer(NHO),
                    PY2=integer(NHO),
                    PX3=integer(NHO),
                    PY3=integer(NHO),
                    PX4=integer(NHO),
                    PY4=integer(NHO))
  house_df<-(data.frame(rec_nr=retvals$POSTS,x1=amin+retvals$PX1/100,y1=bmin+retvals$PY1/100,x2=amin+retvals$PX2/100,y2=bmin+retvals$PY2/100,x3=amin+retvals$PX3/100,y3=bmin+retvals$PY3/100,x4=amin+retvals$PX4/100,y4=bmin+retvals$PY4/100)[1:retvals$NR,])
  
  house<-data.frame(data.table::rbindlist(list(house_df[c("rec_nr","x1","y1")], house_df[c("rec_nr","x2","y2")], house_df[c("rec_nr","x3","y3")],house_df[c("rec_nr","x4","y4")]), use.names=FALSE))
  names(house)<-c("rec_nr","x","y")

#loop through points to construct convex hulls and then lines
for (s in 1:retvals$NR) {
  house_subset <- subset(house,rec_nr==s,select=c(rec_nr,x,y))
  house_subset_cps <- grDevices::chull(house_subset[c("x","y")])
  house_cps <- house[rownames(house_subset[house_subset_cps,]),]
    for (r in 1:4) {
      if (r < 4) {
        r2 = r+1
      } else {
        r2 = 1
      }
      rectangles <- rbind(data.frame(marks=as.factor(s),x0=house_cps[r,]$x,
                                    y0=house_cps[r,]$y,x1=house_cps[r2,]$x,
                                    y1=house_cps[r2,]$y),rectangles)
   }
  }

  
  rectangles_len<-nrow(rectangles)
  rec_rearr<-NULL
  for (rec in 1:rectangles_len) {
    rec_row<-rectangles[rec,]
    if((rec_row$x0 < rec_row$x1)) {
      rec_rearr<-rbind(rec_rearr,
                       data.frame(x_1=rec_row$x0,y_1=rec_row$y0,x_2=rec_row$x1,y_2=rec_row$y1))
    }
    else {
      rec_rearr<-rbind(rec_rearr,data.frame(x_1=rec_row$x1,y_1=rec_row$y1,x_2=rec_row$x0,y_2=rec_row$y0))
    }
  }
  rec_rearr_red<-unique(rec_rearr)
  rec_rearr_red<-cbind(row_id=rownames(rec_rearr_red),rec_rearr_red)
  
  # Attach unique edges back to rectangles
  rec_add<-NULL
  rec_rearr_red_len<-nrow(rec_rearr_red)
  for (rec_ in 1:rectangles_len) {
    rec_row<-rectangles[rec_,]
    for (s in 1:rec_rearr_red_len) {
      rrr_row<-rec_rearr_red[s,]
      if(((rec_row$x0 == rrr_row$x_1) &
          (rec_row$y0 == rrr_row$y_1) &
          (rec_row$x1 == rrr_row$x_2) &
          (rec_row$y1 == rrr_row$y_2)) |
         ((rec_row$x0 == rrr_row$x_2) &
          (rec_row$y0 == rrr_row$y_2) &
          (rec_row$x1 == rrr_row$x_1) &
          (rec_row$y1 == rrr_row$y_1))) {
        rec_add<-rbind(rec_add,data.frame(rec_nr=rec_row$marks,line_nr=rrr_row$row_id))
        #break
      }
    }
  }
  
  # Prepare list of edges per rectangle
  comp_list<-NULL
  rec_add_len<-length(unique(rec_add$rec_nr))
  for (t in 1:rec_add_len) {
    line_nr_list<-c()
    rec_add_sub<-subset(rec_add,subset = rec_nr == t)
    line_nr_list<-append(rec_add_sub$line_nr,list(line_nr_list))
    rec_add_sublist<-list(list(t),(line_nr_list[1:4]))
    comp_list<-append(comp_list,list(rec_add_sublist))
  }
  
  # Search for rectangles with identical edges
  l_comp<-length(comp_list)
  comp_list_start=1
  comp_list_end=0
  while (comp_list_start != comp_list_end) {
    comp_list_start<-comp_list_end
    for (x in 1:l_comp) {
      for (y in 1:l_comp) {
        if (isFALSE(identical(comp_list[[x]][[1]], comp_list[[y]][[1]])) & 
            (length(intersect(comp_list[[x]][[2]],comp_list[[y]][[2]])) > 0)) {
          comp_list[[x]][[1]]<-unique(c(comp_list[[x]][[1]],comp_list[[y]][[1]]))
          comp_list[[x]][[2]]<-unique(c(comp_list[[x]][[2]],comp_list[[y]][[2]]))
        }
      }
      comp_list_end<-length(unlist(comp_list))
    }
  }
  
  # Order lists of adjoining rectangles and remove duplicates
  comp_list_red<-NULL
  for (z in 1:l_comp) {
    comp_list_red<-(c(comp_list_red,list(sort(unlist(comp_list[[z]][[1]])))))
  }
  comp_list_red<-unique(comp_list_red)
  
  # Attach groups of adjoining rectangles back to the individual rectangles
  rec_reattach <- NULL
  for (g in 1:rectangles_len) {
    rec_row<-rectangles[g,]
    for (h in 1:length(comp_list_red)) {
      comp_list_red_row<-comp_list_red[[h]]
      if (length(intersect(rec_row$marks,comp_list_red_row)) > 0) {
        rec_reattach <- rbind(rec_reattach,
                              data.frame(rec_row,group_nr=h,size=length(comp_list_red[[h]])))
        break
      }
    }
  }
  
  X_points <- ppp(house$x, house$y,  checkdup=FALSE, window=owin(c(amin-10,max(a)+10),c(bmin-10,max(b)+10)))
  X_line <- as.psp(rectangles, window=owin(c(amin-10,max(a)+10),c(bmin-10,max(b)+10)))
  
  x_points<-as.SpatialPoints.ppp(X_points)
  x_lines<-as.SpatialLines.psp(X_line)
  
  x_points_shp<-SpatialPointsDataFrame(x_points,data.frame(id=1:length(x_points)))
  x_lines_shp<-SpatialLinesDataFrame(x_lines,data.frame(id=1:length(x_lines),rect_nr=rectangles$marks,group_nr=rec_reattach$group_nr,group_size=rec_reattach$size))
  
  # Define output method
  if (outputmethod=="Spatial") {
    points <- x_points
    lines <- x_lines
  } else if (outputmethod=="shape") {
    points <- x_points_shp
    lines <- x_lines_shp
  } else if (outputmethod=="spatstat") {
    points <- X_points
    lines <- X_line
  }
  return(list(points = points, lines = lines))
}
