##Archaeological analysis=group
##Path_to_postAAR_Fortran=file
##Point_layer=vector
##Minimum_distance = number 2
##Maximum_distance = number 5
##Tolerance_in_locations = number 0.2
##postAAR_points=output vector
##postAAR_lines=output vector

library("spatstat")
library("maptools")
library("sp")
library("data.table")
options(digits=10)

output_method<-"QGIS"

ppp=as(as(Point_layer, "SpatialPoints"),"ppp")
ppp_coords<-coords(ppp)

a<-ppp_coords$x
b<-ppp_coords$y

x<-Minimum_distance
y<-Maximum_distance
z<-Tolerance_in_locations

if (!is.loaded('postaar')) {
dyn.load(Path_to_postAAR_Fortran)
}
NH=0
NR=0
NRE=0

amin<-min(a)
bmin<-min(b)
a<-a-amin
b<-b-bmin

a_<-as.integer(round(a*100))
b_<-as.integer(round(b*100))
z_<-as.integer(round(z*100))
x_<-as.integer(round(x*100))
y_<-as.integer(round(y*100))
NHO=length(a_)
NSO=1000000
retvals<-.Fortran("postaar",
SMIN=as.integer(x_),
SMAX=as.integer(y_),
TOL=as.integer(z_),
NH=as.integer(NH),
NR=as.integer(NR),
NHO=as.integer(NHO),
NRE=as.integer(NRE),
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

house<-data.frame(rbindlist(list(house_df[c("rec_nr","x1","y1")], house_df[c("rec_nr","x2","y2")], house_df[c("rec_nr","x3","y3")],house_df[c("rec_nr","x4","y4")])))
names(house)<-c("rec_nr","x","y")

rectangles <- NULL

for (s in 1:retvals$NR) {
house_subset<-subset(house,rec_nr==s,select=c(rec_nr,x,y))
house_subset_cps <- chull(house_subset[c("x","y")])
house_cps<-house[rownames(house_subset[house_subset_cps,]),]
for (r in 1:4) {
if (r < 4) {
r2=r+1
} else {
r2=1
}
rectangles <- rbind(data.frame(marks=as.factor(s),x0=house_cps[r,]$x,
y0=house_cps[r,]$y,x1=house_cps[r2,]$x,
y1=house_cps[r2,]$y),rectangles)
}
}

# rearrange order of vertices of rectangle and reduce edges to unique ones
pa<-rectangles
pa_len<-nrow(pa)
rec_rearr<-NULL
for (rec in 1:pa_len) {
  pa_row<-pa[rec,]
  if((pa_row$x0 < pa_row$x1)) {
    rec_rearr<-rbind(rec_rearr,
                     data.frame(x_1=pa_row$x0,y_1=pa_row$y0,x_2=pa_row$x1,y_2=pa_row$y1))
  }
  else {
    rec_rearr<-rbind(rec_rearr,data.frame(x_1=pa_row$x1,y_1=pa_row$y1,x_2=pa_row$x0,y_2=pa_row$y0))
  }
}
rec_rearr_red<-unique(rec_rearr)
rec_rearr_red<-cbind(row_id=rownames(rec_rearr_red),rec_rearr_red)

# Attach unique edges back to rectangles
pa_add<-NULL
rec_rearr_red_len<-nrow(rec_rearr_red)
for (rec_ in 1:pa_len) {
  pa_row<-pa[rec_,]
  for (s in 1:rec_rearr_red_len) {
    rrr_row<-rec_rearr_red[s,]
  if(((pa_row$x0 == rrr_row$x_1) &
      (pa_row$y0 == rrr_row$y_1) &
      (pa_row$x1 == rrr_row$x_2) &
      (pa_row$y1 == rrr_row$y_2)) |
      ((pa_row$x0 == rrr_row$x_2) &
      (pa_row$y0 == rrr_row$y_2) &
      (pa_row$x1 == rrr_row$x_1) &
      (pa_row$y1 == rrr_row$y_1))) {
    pa_add<-rbind(pa_add,data.frame(rec_nr=pa_row$marks,line_nr=rrr_row$row_id))
    break
  }
}
}

# Prepare list of edges per rectangle
comp_list<-NULL
pa_add_len<-length(unique(pa_add$rec_nr))
for (t in 1:pa_add_len) {
  line_nr_list<-c()
  pa_add_sub<-subset(pa_add,subset = rec_nr == t)
  line_nr_list<-append(pa_add_sub$line_nr,list(line_nr_list))
  pa_add_sublist<-list(list(t),(line_nr_list[1:4]))
  comp_list<-append(comp_list,list(pa_add_sublist))
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
pa_reattach <- NULL
for (g in 1:pa_len) {
  pa_row<-pa[g,]
  for (h in 1:length(comp_list_red)) {
    comp_list_red_row<-comp_list_red[[h]]
    if (length(intersect(pa_row$marks,comp_list_red_row)) > 0) {
      pa_reattach <- rbind(pa_reattach,
                           data.frame(pa_row,group_nr=h,size=length(comp_list_red[[h]])))
    break
    }
  }
}

X_points <- ppp(house$x, house$y,  checkdup=FALSE, window=owin(c(amin+min(a_/100)-10,amin+max(a_/100)+10),c(bmin+min(b_/100)-10,max(bmin+b_/100)+10)))
X_line <- as.psp(rectangles, window=owin(c(amin+min(a_/100)-10,amin+max(a_/100)+10),c(bmin+min(b_/100)-10,max(bmin+b_/100)+10)))

x_points<-as.SpatialPoints.ppp(X_points)
x_lines<-as.SpatialLines.psp(X_line)

x_points_shp<-SpatialPointsDataFrame(x_points,data.frame(id=1:length(x_points)))
x_lines_shp<-SpatialLinesDataFrame(x_lines,data.frame(id=1:length(x_lines),rect_nr=rectangles$marks,group_nr=pa_reattach$group_nr,group_size=pa_reattach$size))

# Define output method
if (output_method=="QGIS") {
  postAAR_points=x_points_shp
  postAAR_lines=x_lines_shp
} else if (output_method=="Spatial") {
  return(list(points=x_points,lines=x_lines))
} else if (output_method=="shape") {
  return(list(points=x_points_shp,lines=x_lines_shp))
} else if (output_method=="spatstat") {
  return(list(points=X_points,lines=X_line))
}


comp_list_red_len<-lengths(comp_list_red)
comp_list_red_sum<-as.data.frame(table(comp_list_red_len))
colnames(comp_list_red_sum)<-c("rec_num","count")
comp_list_red_output<-NULL
for (v in 1:nrow(comp_list_red_sum)) {
  comp_list_red_output<-paste(sep = "\n", comp_list_red_output,
                          paste(sep = " "," ", comp_list_red_sum$rec_num[v]," ",comp_list_red_sum$count[v]))
}
>cat(sep = '',"\n", retvals$NH, " postholes were read in.\n\nWith the given parameters (minimum distance: ", retvals$SMIN/100, ", maximum distance: ", retvals$SMAX/100, ", tolerance of location: ", retvals$TOL/100, "), postAAR found ", retvals$NR, " rectangles.\nThe expected number was ", retvals$NRE,".\n\nNumber and count of joined rectangles:",comp_list_red_output)
