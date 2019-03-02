

##input similarity matrix of community, geodistance matrix, geodistance scale##


disdecay<-function(commsimi,geodis, geoscale=c(0,Inf),nperm=999){
  Y<-commsimi[row(commsimi)>col(commsimi)]
  X<-geodis[row(geodis)>col(geodis)]
  minY<-min(Y[Y>0])
  Y[Y<=0]<-minY
  logY<-log10(Y)
  logX<-log10(X+1)
  logY1<-logY[X>=geoscale[1] & X<geoscale[2]]
  logX1<-logX[X>=geoscale[1] & X<geoscale[2]]
  
  logX1<-cbind(rep(1,length(logX1)),logX1)
  
  XX <- crossprod(logX1)
  XX <- solve(XX)
  
  # will need to calculate Xy for each permutation
  XY <- crossprod(logX1, logY1)

   # regression coefficients
  b <- XX %*% XY
  #slop
  slop.obs<-b[2]
  
  perm.slope<-sapply(1:nperm,function(i){
    newY <- ecodist::full(logY)
    newSample <- sample(nrow(newY))
    newY <- newY[newSample, newSample]
    logYnew <- ecodist::lower(newY)
    logYnew<-logYnew[X>=geoscale[1] & X<geoscale[2]]
    
    
    # will need to calculate Xy for each permutation
    XY <- crossprod(logX1, logYnew)
    
    # regression coefficients
    b <- XX %*% XY
    #slop
    b[2]
  })
  
  permSlope.sd<-sd(perm.slope)
  slope.rank<-rank(c(slop.obs,perm.slope),ties.method="first")[1]
  ttest.pvalue<-t.test(perm.slope,mu=slop.obs)$p.value
  
  c(slope.obs=slop.obs,permSlp.mean=mean(perm.slope),permSlope.sd=permSlope.sd,slope.rank=slope.rank,nperm=nperm,ttest.pvalue=ttest.pvalue)
}

disdecay2<-function(commsimi,geodis, geoscale=c(0,Inf),nperm=999){
  Y<-commsimi[row(commsimi)>col(commsimi)]
  X<-geodis[row(geodis)>col(geodis)]
  minY<-min(Y[Y>0])
  Y[Y<=0]<-minY
  logY<-log10(Y)
  logX<-log10(X+1)
  logY1<-logY[X>=geoscale[1] & X<geoscale[2]]
  logX1<-logX[X>=geoscale[1] & X<geoscale[2]]
  
  logX1<-cbind(rep(1,length(logX1)),logX1)
  
  XX <- crossprod(logX1)
  XX <- solve(XX)
  
  # will need to calculate Xy for each permutation
  XY <- crossprod(logX1, logY1)
  
  # regression coefficients
  b <- XX %*% XY
  #slop
  slop.obs<-b[2]
  
  perm.slope<-sapply(1:nperm,function(i){
    newY <- ecodist::full(logY)
    newSample <- sample(nrow(newY))
    newY <- newY[newSample, newSample]
    logYnew <- ecodist::lower(newY)
    logYnew<-logYnew[X>=geoscale[1] & X<geoscale[2]]
    
    
    # will need to calculate Xy for each permutation
    XY <- crossprod(logX1, logYnew)
    
    # regression coefficients
    b <- XX %*% XY
    #slop
    b[2]
  })
  result<-c(slop.obs,perm.slope)
  names(result)<-c("observed",sprintf("perm%s",1:nperm))
  result
}


#test

canberradis<-read.csv("/Users/linwei/Dropbox/global water2/6-distance decay/canberra distance.csv",row.names = 1,header = T)

load("~/Dropbox/global water2/4-DCA pcoa/GW16S.taxadist.rda")
load("~/Dropbox/global water2/4-DCA pcoa/GW16S.unifrac.rda")

sorensen=as.matrix(betaTD$sorensen);bray=as.matrix(betaTD$bray);
unfd_1=as.matrix(ufn$d_1);unfUW=as.matrix(ufn$d_UW)

sorensen=sorensen[match(row.names(canberradis),row.names(sorensen)),match(row.names(canberradis),row.names(sorensen))]
bray=bray[match(row.names(canberradis),row.names(bray)),match(row.names(canberradis),row.names(bray))]
unfd_1=unfd_1[match(row.names(canberradis),row.names(unfd_1)),match(row.names(canberradis),row.names(unfd_1))]
unfUW=unfUW[match(row.names(canberradis),row.names(unfUW)),match(row.names(canberradis),row.names(unfUW))]


allcommdist<-list(sorensen=sorensen,bray=bray,
                 unfd_1=unfd_1,unfUW=unfUW,canberradis=canberradis)




geodis<-read.csv("~/Dropbox/global water2/6-distance decay/geodistance.correct.csv",row.names = 1,header = T)

geodis<-geodis[match(row.names(canberradis),row.names(geodis)),match(row.names(canberradis),row.names(geodis))]

##within city level, different similarity index##

withincity.result<-sapply(1:5,function(i){
  commdis<-as.matrix(allcommdist[[i]])
  commdis<-commdis[match(row.names(geodis),row.names(commdis)),match(row.names(geodis),row.names(commdis))]
  commsimilar<-1-commdis
  disdecay(commsimi=commsimilar,geodis=geodis,geoscale=c(0,1E5))
})

colnames(withincity.result)<-names(allcommdist)

write.csv(withincity.result,"within site 1E5 distance decay_matrix permed.csv")

##across city, within continent##

acrosscity.result<-sapply(1:5,function(i){
  commdis<-as.matrix(allcommdist[[i]])
  commdis<-commdis[match(row.names(geodis),row.names(commdis)),match(row.names(geodis),row.names(commdis))]
  commsimilar<-1-commdis
  disdecay(commsimi=commsimilar,geodis=geodis,geoscale=c(1E5,5E6))
})

colnames(acrosscity.result)<-names(allcommdist)

write.csv(acrosscity.result,"across site 5E6 distance decay_permued matrix.csv")


##across continent##

acrossconti.result<-sapply(1:5,function(i){
  commdis<-as.matrix(allcommdist[[i]])
  commdis<-commdis[match(row.names(geodis),row.names(commdis)),match(row.names(geodis),row.names(commdis))]
  commsimilar<-1-commdis
  disdecay(commsimi=commsimilar,geodis=geodis,geoscale=c(5E6,Inf))
})

colnames(acrossconti.result)<-names(allcommdist)

write.csv(acrossconti.result,"across continent distance decay_permued matrix.csv")

##all
all.result<-sapply(1:5,function(i){
  commdis<-as.matrix(allcommdist[[i]])
  commdis<-commdis[match(row.names(geodis),row.names(commdis)),match(row.names(geodis),row.names(commdis))]
  commsimilar<-1-commdis
  disdecay(commsimi=commsimilar,geodis=geodis,geoscale=c(0,Inf))
})

colnames(all.result)<-names(allcommdist)

write.csv(all.result,"global distance decay_permued matrix.csv")

#within city level, permutation result

withincity.result2<-sapply(1:5,function(i){
  commdis<-as.matrix(allcommdist[[i]])
  commdis<-commdis[match(row.names(geodis),row.names(commdis)),match(row.names(geodis),row.names(commdis))]
  commsimilar<-1-commdis
  disdecay2(commsimi=commsimilar,geodis=geodis,geoscale=c(0,1E5))
})
colnames(withincity.result2)<-names(allcommdist)
write.csv(withincity.result2,"within site slope perms_by matrix.csv")

#across city level, permutation result

acrosscity.result2<-sapply(1:5,function(i){
  commdis<-as.matrix(allcommdist[[i]])
  commdis<-commdis[match(row.names(geodis),row.names(commdis)),match(row.names(geodis),row.names(commdis))]
  commsimilar<-1-commdis
  disdecay2(commsimi=commsimilar,geodis=geodis,geoscale=c(1E5,5E6))
})
colnames(acrosscity.result2)<-names(allcommdist)
write.csv(acrosscity.result2,"across site slope perms_by matrix.csv")

#across continent level, permutation result

acrossconti.result2<-sapply(1:5,function(i){
  commdis<-as.matrix(allcommdist[[i]])
  commdis<-commdis[match(row.names(geodis),row.names(commdis)),match(row.names(geodis),row.names(commdis))]
  commsimilar<-1-commdis
  disdecay2(commsimi=commsimilar,geodis=geodis,geoscale=c(5E6,Inf))
})
colnames(acrossconti.result2)<-names(allcommdist)
write.csv(acrossconti.result2,"across conti slope perms_by matrix.csv")

##all
all.result2<-sapply(1:5,function(i){
  commdis<-as.matrix(allcommdist[[i]])
  commdis<-commdis[match(row.names(geodis),row.names(commdis)),match(row.names(geodis),row.names(commdis))]
  commsimilar<-1-commdis
  disdecay2(commsimi=commsimilar,geodis=geodis,geoscale=c(0,Inf))
})
colnames(all.result2)<-names(allcommdist)
write.csv(all.result2,"global slope perms_by matrix.csv")

##get rank of slope difference#compare slope in each scale with overall slope#
##within city v.s. overall slope

diffslope<-withincity.result2-all.result2
rank.within<-sapply(1:5,function(i){
  obsperm<-diffslope[,i]
  rank(obsperm,ties.method="first")[1]
})

diffslope2<-acrosscity.result2-all.result2
rank.acrosscity<-sapply(1:5,function(i){
  obsperm<-diffslope2[,i]
  rank(obsperm,ties.method="first")[1]
})

diffslope3<-acrossconti.result2-all.result2
rank.acrossconti<-sapply(1:5,function(i){
  obsperm<-diffslope3[,i]
  rank(obsperm,ties.method="first")[1]
})

rank.result<-cbind(rank.within,rank.acrosscity,rank.acrossconti)
row.names(rank.result)<-names(allcommdist)
write.csv(rank.result,"compare slope with overall.rank_by matrix.csv")



