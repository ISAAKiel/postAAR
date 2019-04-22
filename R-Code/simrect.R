simrect<-function(a,b,c,d,e,f,g,rp,h,k) {
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
    cr<-sample(2:h, size = 1)
    rows<-sample(2:k, size = 1)
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
  
  cat(sep = '',"\n",c, " houses with ", rectangle_nr, " rectangles were simulated.\n\n")

if (rp>0) {
  # generate additional randompoints
  rand_x<-runif(rp, a, (b+2*h*e)) # x-coordinates
  rand_y<-runif(rp, a, (b+2*h*e)) # y-coordinates
  
  houses<-rbind(data.frame(rec_nr=NA,x=rand_x,y=rand_y),houses)
  }
  
  # create an object of type pointpattern in spatstat
  X <- ppp(houses$x, houses$y, c(a,(b+2*h*e)),c(a,(b+2*h*e)))

  # create an object of type linepattern in spatstat
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
  
  list(points=X_jittered,lines=X_line)

}
