# rectangles aus postAAR = pa
# rearrange order of vertices of rectangle and reduce edges to unique ones
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
head(pa_reattach)
comp_list_red_len<-lengths(comp_list_red)
str(comp_list_red_len)
comp_list_red_sum<-as.data.frame(table(comp_list_red_len))
colnames(comp_list_red_sum)<-c("rec_num","count")
comp_list_red_sum
comp_list_red_output<-NULL
for (v in 1:nrow(comp_list_red_sum)) {
  comp_list_red_output<-paste(sep = "\n", comp_list_red_output,
                          paste(sep = "","number of rectangles: ",comp_list_red_sum$rec_num[v], ", count: ",comp_list_red_sum$count[v]))
}
cat(comp_list_red_output)
