##Archaeological analysis=group
##Lower_left_corner = number 0
##Upper_right_corner = number 200
##Number_of_houses = number 50
##Maximum_number_of_rows_and_crossrows = number 5
##Minimum_distance = number 2
##Maximum_distance = number 5
##Point_jitter = number 0.2
##Preservation_percentage = number 100
##Additional_random_points = number 1000
##simpostAAR_points=output vector
##simpostAAR_lines=output vector

library("spatstat")
library("maptools")
library("sp")
library("data.table")

a<-Lower_left_corner
b <-Upper_right_corner
c<-Number_of_houses
d<-Minimum_distance
e<-Maximum_distance
f<-Point_jitter
g<-Preservation_percentage
h<-Maximum_number_of_rows_and_crossrows
rp<-Additional_random_points
houses<-NULL
house_sim<-NULL
rectangle_nr<-0
for (i in 1:c) {
house<-NULL
rec_nr<-i
x<-sample((a+h*e):(b+h*e), size = 1)
y<-sample((a+h*e):(b+h*e), size = 1)
w<-sample(d:e, size = 1)
l<-sample(d:e, size = 1)

if (h<3) {
cr=2
rows=2
} else {
cr<-sample(2:h, size = 1)
rows<-sample(2:h, size = 1)
}
theta<-sample(1:180, size =1)
rectangle_nr<-(cr-1)*(rows-1)+rectangle_nr

for (j in 1:cr) {
for (m in 1:rows) {
tempx<- -(-0.5*(rows+1)+m)*w
tempy<- -(-0.5*(cr+1)+j)*l
rotatedx <- tempx*cos(theta) - tempy*sin(theta);
rotatedy <- tempx*sin(theta) + tempy*cos(theta);
ax_rot<- rotatedx + x
ay_rot<- rotatedy + y
house<-rbind(data.frame(rec_nr,x=ax_rot,y=ay_rot),house)
}
}
houses<-rbind(house,houses)
house_cps <- chull(house$x,house$y)
cps_nr <-length(house_cps)
for (r in 1:cps_nr) {
if (r < cps_nr) {
r2=r+1
} else {
r2=1
}
house_sim <- rbind(data.frame(marks=as.factor(rec_nr),x0=house[house_cps[r],]$x,
y0=house[house_cps[r],]$y,x1=house[house_cps[r2],]$x,
y1=house[house_cps[r2],]$y),house_sim)
}
}
if (rp>0) {
# generate additional randompoints
rand_x<-runif(rp, a, (b+2*h*e)) # x-coordinates
rand_y<-runif(rp, a, (b+2*h*e)) # y-coordinates

houses<-rbind(data.frame(rec_nr=NA,x=rand_x,y=rand_y),houses)
}
X <- ppp(houses$x, houses$y, c(a,(b+2*h*e)),c(a,(b+2*h*e)))
X_line <- as.psp(house_sim, window=owin(c(a,(b+2*h*e)),c(a,(b+2*h*e))))
if (g<100) {
# thin the dataset
X_thinned <- rthin(X,g/100)
} else {
X_thinned<-X
}

if (f>0) {
# jitter the points
X_jittered <-rjitter(X_thinned,f)
} else {
X_jittered<-X_thinned
}

simrect_list<-list(points=X_jittered,lines=X_line)

sim_points<-as.SpatialPoints.ppp(simrect_list$points)
simpoints_shp<-SpatialPointsDataFrame(sim_points,data.frame(id=1:length(sim_points)))

sim_lines<-as.SpatialLines.psp(simrect_list$lines)
simlines_shp<-SpatialLinesDataFrame(sim_lines,data.frame(id=1:length(sim_lines)))
simpostAAR_points=simpoints_shp
simpostAAR_lines=simlines_shp
>cat(sep = '',"\n",c, " houses with ", rectangle_nr, " rectangles were simulated.\n\n")
