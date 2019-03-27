postAAR<- function(a,b,x,y,z) {
  if (!is.loaded('postaar')) {
    dyn.load("postAAR_subroutine.so")
  }
  NH=0
  NR=0
  NRE=0
  NHO=length(a)
  NSO=1000000
  retvals<-.Fortran("postaar",
                    SMIN=as.integer(x*100),
                    SMAX=as.integer(y*100),
                    TOL=as.integer(z*100),
                    NH=as.integer(NH),
                    NR=as.integer(NR),
                    NHO=as.integer(NHO),
                    NRE=as.integer(NRE),
                    NSO=as.integer(NSO),
                    IX=as.integer(a*100),
                    IY=as.integer(b*100),
                    POSTS=integer(NHO),
                    PX1=integer(NHO),
                    PY1=integer(NHO),
                    PX2=integer(NHO),
                    PY2=integer(NHO),
                    PX3=integer(NHO),
                    PY3=integer(NHO),
                    PX4=integer(NHO),
                    PY4=integer(NHO))
  cat(sep = '',"\n", retvals$NH, " postholes were read in.\n\nWith the given parameters (minimal distance: ", retvals$SMIN/100, ", maximal distance: ", retvals$SMAX/100, ", tolerance of angle: ", retvals$TOL/100, "), postAAR found ", retvals$NR, " rectangles.\nThe expected number was ", retvals$NRE,".\n\n")
  house_df<-(data.frame(rec_nr=retvals$POSTS,x1=retvals$PX1/100,y1=retvals$PY1/100,x2=retvals$PX2/100,y2=retvals$PY2/100,x3=retvals$PX3/100,y3=retvals$PY3/100,x4=retvals$PX4/100,y4=retvals$PY4/100)[1:retvals$NR,])

  house<-data.frame(rbindlist(list(house_df[c("rec_nr","x1","y1")], house_df[c("rec_nr","x2","y2")], house_df[c("rec_nr","x3","y3")],house_df[c("rec_nr","x4","y4")])))
  names(house)<-c("rec_nr","x","y")
  
    # create an object of type pointpattern in spatstat
  X_points <- ppp(house$x, house$y,  window=owin(c(min(a)-10,max(a)+10),c(min(b)-10,max(b)+10)))

rectangles <- NULL
#loop through points to construct convex hulls and then lines
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
  
  # create an object of type linepattern in spatstat
  X_line <- as.psp(rectangles, window=owin(c(min(a)-10,max(a)+10),c(min(b)-10,max(b)+10)))
  
  return(list(points=X_points,lines=X_line))
}