distLLE<-function(latitude,longitude,elevation=NA,site.name=NA,short=TRUE)
{
# To calculate distance of two points according to longitude, latitude and elevation.
# by Daliang Ning (ningdaliang@gmail.com) on 2015.2.17
# if short=TRUE, the points are close to each other,dist2=H2+d2.
  library(geosphere)
  num=length(latitude)
  if(is.na(site.name[1])){site.name=paste("S",1:num,sep="")}
  x=data.frame(longitude,latitude)
  d=distm(x)
  if(is.na(elevation[1]))
  {
    dist=d
  }else{
    H=as.matrix(dist(elevation))
    if(short)
    {
      dist=(d^2+H^2)^0.5
    }else{
      h=matrix(pmin(elevation[row(matrix(,num,num))],elevation[col(matrix(,num,num))]),nrow=num)
      R=6378137
      dist=(((2*(R+h)*sin(0.5*d/R)*sin(0.5*d/R)+H)^2)+(((R+h)*sin(d/R))^2))^0.5
    }
  }
  rownames(dist)=site.name
  colnames(dist)=site.name
  dist
}

