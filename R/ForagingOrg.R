best.r<-function(x,y){
  if (length(x)==0 | length(y)==0){
    return(NA)
  }else{
  order<-c(1:length(x))
  corx<-abs(cor(x,order))
  cory<-abs(cor(y,order))
  return(max(corx,cory))
  }
}

meanITD<-function(x,y){
  itd_matrix<-dist(data.frame(lng=x,lat=y),diag=TRUE,upper=TRUE)
  itd_matrix2<-as.matrix(itd_matrix,nrows=sqrt(length(itd_matrix)))
  itd_matrix3<-itd_matrix2[-nrow(itd_matrix2),-1]
  path_length<-sum(diag(itd_matrix3))
  if (nrow(itd_matrix2)>2){
    mean_itd<-mean(diag(itd_matrix3))
  }
  if (nrow(itd_matrix2)==2){
    mean_itd<-as.numeric(itd_matrix3)
  }
  return(mean_itd)
}

PAO<-function(x,y){
  if (length(x)>3){
    opt_itd<-numeric()
    itd_matrix<-dist(data.frame(lng=x,lat=y),diag=TRUE,upper=TRUE)
    itd_matrix2<-as.matrix(itd_matrix,nrows=sqrt(length(itd_matrix)))
    itd_tsp<-TSP::TSP(itd_matrix)
    tour<-PairViz::order_tsp(itd_matrix,cycle=FALSE)
    for (ot in 1:(length(tour)-1)){
      opt_itd[ot]<-itd_matrix2[tour[ot],tour[ot+1]]
    }
    opt_tour_length<-sum(opt_itd)
    itd_matrix3<-itd_matrix2[-nrow(itd_matrix2),-1]
    scanpath_length<-sum(diag(itd_matrix3))
    PAO<-(scanpath_length / opt_tour_length -1)*100
  }else{PAO<-NA}
  return(PAO)
}

number.intersections<-function(x,y){
  intersection_counter<-0
  if(length(x)>3){
    for (p in 1:(length(x)-3)){
      if (x[p]==x[p+1]){x[p+1]<-x[p+1]+0.01}
      for (p2 in (p+2):(length(x)-1)){

        if (x[p2]==x[p2+1]){x[p2+1]<-x[p2+1]+0.01}

        x1a<-min(x[p],x[p+1])
        x2a<-max(x[p],x[p+1])
        x1b<-min(x[p2],x[p2+1])
        x2b<-max(x[p2],x[p2+1])

        if (x1a==x[p]){
          y1a<-y[p]
          y2a<-y[p+1]
        }
        if (x1a==x[p+1]){
          y1a<-y[p+1]
          y2a<-y[p]
        }
        if (x1b==x[p2]){
          y1b<-y[p2]
          y2b<-y[p2+1]
        }
        if (x1b==x[p2+1]){
          y1b<-y[p2+1]
          y2b<-y[p2]
        }
        xi<-max(x1a,x1b)
        xf<-min(x2a,x2b)

        r1<-approxfun(x=c(x1a,x2a),y=c(y1a,y2a))
        r2<-approxfun(x=c(x1b,x2b),y=c(y1b,y2b))

        if(is.na(r1(xi))==FALSE & is.na(r2(xi))==FALSE & is.na(r1(xf))==FALSE & is.na(r2(xf))==FALSE & (sign(r1(xi)-r2(xi)) != sign(r1(xf)-r2(xf)))){
          intersection_counter<-intersection_counter+1
        }
      }
    }
  }
  intersection_rate<-intersection_counter/length(x)
  intersection_list<-list(intersection_counter,intersection_rate)
  names(intersection_list)<-c("Number.of.intersections","Intersection.rate")
  return(intersection_list)
}

Q.score<-function(n.collected.targets,n.total.targets,rt){
  Q<-n.collected.targets^2/(n.total.targets*rt)
  return(Q)
}

