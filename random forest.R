
env.file="meta_plant level.csv"
com.file="OTU table_plant level.csv"


# internal use. not published yet. please at least acknowledge if you used.
#install.packages("ieggr",repos = c("ftp://gr:gr!123@129.15.40.254","http://cran.us.r-project.org"))
# In some cases, you need to install it in R rather than RStudio.
library(ieggr)
library(parallel)
n<-detectCores()

comm=t(lazyopen(com.file))

env.sample<-lazyopen(env.file)
head(env.sample)
group<-env.sample$Continent


samp.all=rownames(env.sample)
check.samp=match.name(name.check = samp.all,rn.list = list(comm=comm))
comm<-check.samp$comm;sum(row.names(comm)==row.names(env.sample))

##random forest all data##
test1<-colSums(comm>0)
Group<-factor(group,levels=unique(group))
test2<-sapply(split(as.data.frame(comm),Group),colSums)
Test2<-rowSums(test2>0)  ##detected in how many continents

otuid<-which(test1>=55 & Test2>=6) ##detected in 20% samples, and in every continent## 
commused<-comm[,otuid]

env.sample$Influent.BOD..mg.L.[env.sample$Influent.BOD..mg.L.>800]<-NA   ##remove outliers in response variable
yenv<-env.sample$Influent.BOD..mg.L.

sampleid<-which(!is.na(yenv))
yenvused<-yenv[sampleid]
xcommused<-commused[sampleid,]

library(randomForest)
randomf<-randomForest(x=xcommused,y=yenvused,ntree=500,proximity=T,importance = T)
print(randomf)
a<-randomf$predicted
crosspredict<-data.frame(observedBOD=yenvused,crosspredictedBOD=a)
plot(yenvused,a)
save.file(crosspredict,filename = "cross predict_BOD.otu",folder = save.wd)

b<-randomf$importance
save.file(b,filename = "importance_BOD.otu",folder = save.wd)

c<-randomf$proximity
save.file(c,filename = "proximity_BOD.otu",folder = save.wd)

##psedo r2, to indicate prediction strength#
ymean<-mean(crosspredict$observedBOD);
sum((crosspredict$observedBOD-crosspredict$crosspredictedBOD)^2)/sum((crosspredict$observedBOD-ymean)^2)-1


##sludge age
env.sample$Sludge.Age..days.[env.sample$Sludge.Age..days.>35]<-NA
yenv<-env.sample$Sludge.Age..days.
sampleid<-which(!is.na(yenv))
yenvused<-yenv[sampleid]
xcommused<-commused[sampleid,]

randomf<-randomForest(x=xcommused,y=yenvused,ntree=500,proximity=T,importance = T)
print(randomf)
a<-randomf$predicted
crosspredict<-data.frame(observedSage=yenvused,crosspredictedSage=a)
plot(yenvused,a)
save.file(crosspredict,filename = "cross predict_Sludge age.otu",folder = save.wd)

b<-randomf$importance
save.file(b,filename = "importance_Sludge age.otu",folder = save.wd)

c<-randomf$proximity
save.file(c,filename = "proximity_Sludge age.otu",folder = save.wd)


##psedo r2#
ymean<-mean(crosspredict$observedSage);

sum((crosspredict$observedSage-crosspredict$crosspredictedSage)^2)/sum((crosspredict$observedSage-ymean)^2)-1


##Temperature
yenv<-env.sample$Air.temperature.Mean.annual
sampleid<-which(!is.na(yenv))
yenvused<-yenv[sampleid]
xcommused<-commused[sampleid,]

trainid<-which(env.sample$Continent %in% c("North America","South America"))
testid<-which(env.sample$Continent=="Asia")

xcommused<-commused[trainid,];yenvused<-yenv[trainid]
randomf<-randomForest(x=xcommused,y=yenvused,ntree=500,proximity=T,importance = T)
print(randomf)
a<-randomf$predicted
crosspredict<-data.frame(observedTemp=yenvused,crosspredictedTemp=a)
plot(yenvused,a)
save.file(crosspredict,filename = "cross predict_Temp.otu",folder = save.wd)

b<-randomf$importance
save.file(b,filename = "importance_Temp.otu",folder = save.wd)

c<-randomf$proximity
save.file(c,filename = "proximity_Temp.otu",folder = save.wd)

trainid<-which(env.sample$Continent %in% c("North America","South America"))
testid<-which(env.sample$Continent=="Asia")

xcommused<-commused[trainid,];yenvused<-yenv[trainid]
randomf<-randomForest(x=xcommused,y=yenvused,ntree=500,proximity=T,importance = T)

testcommused<-commused[testid,]
testTemp.used<-yenv[testid]
d<-predict(randomf,testcommused)
testpredict<-data.frame(observedTemp=testTemp.used,testpredicted.temp=d)
save.file(testpredict,filename = "test predict_Temp",folder = save.wd)
cor(d,testTemp.used)
plot(testTemp.used,d)

##plot
library(ggplot2)
dat2<-data.frame(x=c(2,29),y=c(2,29))
testpredict$testpredicted.temp
ggplot(testpredict,aes(x=observedTemp,y=testpredicted.temp))+
  geom_point(alpha=0.6,colour="blue")+
  geom_segment(aes(x = 5, y = 5, xend = 28, yend = 28),colour="red")+
  xlab("Observed temperature") +
  ylab("Predicted temperature") +
  xlim(4,28)+
  ylim(4,28)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

##psedo r2#
ymean<-mean(testpredict$observedTemp);

sum((testpredict$observedTemp-testpredict$testpredicted.temp)^2)/sum((testpredict$observedTemp-ymean)^2)-1

##temp
crosspredict$observedTemp; crosspredict$crosspredictedTemp
ggplot(crosspredict,aes(x=observedTemp,y=crosspredictedTemp))+
  geom_point(alpha=0.6,colour="blue")+
  geom_segment(aes(x = 4, y = 4, xend = 28, yend = 28),colour="red")+
  xlab("Observed temperature") +
  ylab("Predicted temperature") +
  xlim(2,30)+
  ylim(2,30)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

##psedo r2#
ymean<-mean(crosspredict$observedTemp);

sum((crosspredict$observedTemp-crosspredict$crosspredictedTemp)^2)/sum((crosspredict$observedTemp-ymean)^2)-1


##psedo r2# y is observed, x is predicted
ymean<-mean(y)
sum((y-x)^2)/sum((y-ymean)^2)-1


###########exclude geographical neighbors###############################
source("distLLE.r")
geodis=distLLE(latitude = env.sample$Latitude.Google,longitude = env.sample$Longitude.Google,site.name = rownames(env.sample),short = FALSE)


geoexclude<-function(xcomm,yenv,geodis,distance) {
  out<-t(parSapply(c1,1:nrow(xcomm), function (i) {
    remain<-which(geodis[i,]>=distance)
    names(remain)<-NULL
    remain2<-sample(c(1:nrow(xcomm))[-i],length(remain))
    randomf<-randomForest(x=xcomm[remain,],y=yenv[remain],ntree=500)
    excpredicted<-predict(randomf,newdata=xcomm[i,])
    randomf2<-randomForest(x=xcomm[remain2,],y=yenv[remain2],ntree=500)
    randompredicted<-predict(randomf2,newdata=xcomm[i,])
    observed<-yenv[i]
    c(observed,excpredicted,randompredicted)
  }))
  row.names(out)<-row.names(xcomm)
  colnames(out)<-c("observed","exclude.predict","rdexclu.predict")
  out
}



##sludge age genus##
env.sample$Sludge.Age..days.[env.sample$Sludge.Age..days.>35]<-NA
yenv<-env.sample$Sludge.Age..days.
sampleid<-which(!is.na(yenv))
yenvused<-yenv[sampleid]
xcommused<-commused[sampleid,]
geodisused<-geodis[sampleid,sampleid]
threshold<-c(10,100,500,1000,10000,50000,100000,1000000,5000000)

no_cores=16
c1<-makeCluster(no_cores)
clusterExport(c1,varlist=c("geodisused","threshold","xcommused","yenvused"))
clusterEvalQ(c1, library(randomForest))
diss<-sapply(1:length(threshold),function(i){
  message("Now i=",i," in ",length(threshold),". ",date())
  dis<-threshold[i]
  geoexclude(xcomm=xcommused,yenv=yenvused,geodis=geodisused,distance=dis)
})
colnames(diss)<-c("10m","100m","500m","1000m","10000m","50000m","100000m","1000000m","5000000m")
save.file(diss,filename = "Plant level_Sludge genus level_exclude",folder = save.wd)
stopCluster(c1)

##BOD genus level##

env.sample$Influent.BOD..mg.L.[env.sample$Influent.BOD..mg.L.>800]<-NA
yenv<-env.sample$Influent.BOD..mg.L.

sampleid<-which(!is.na(yenv))
yenvused<-yenv[sampleid]
xcommused<-commused[sampleid,]
geodisused<-geodis[sampleid,sampleid]
threshold<-c(10,100,500,1000,10000,50000,100000,1000000,5000000)

no_cores=16
c1<-makeCluster(no_cores)
clusterExport(c1,varlist=c("geodisused","threshold","xcommused","yenvused"))
clusterEvalQ(c1, library(randomForest))
diss<-sapply(1:length(threshold),function(i){
  message("Now i=",i," in ",length(threshold),". ",date())
  dis<-threshold[i]
  geoexclude(xcomm=xcommused,yenv=yenvused,geodis=geodisused,distance=dis)
})
colnames(diss)<-c("10m","100m","500m","1000m","10000m","50000m","100000m","1000000m","5000000m")
save.file(diss,filename = "plant level_bod genus level_exclude",folder = save.wd)
stopCluster(c1)

## temperature##
yenv<-env.sample$Air.temperature.Mean.annual
sampleid<-which(!is.na(yenv))
yenvused<-yenv[sampleid]
xcommused<-commused[sampleid,]
geodisused<-geodis[sampleid,sampleid]
threshold<-c(10,100,500,1000,10000,50000,100000,1000000,5000000)

no_cores=16
c1<-makeCluster(no_cores)
clusterExport(c1,varlist=c("geodisused","threshold","xcommused","yenvused"))
clusterEvalQ(c1, library(randomForest))
diss<-sapply(1:length(threshold),function(i){
  message("Now i=",i," in ",length(threshold),". ",date())
  dis<-threshold[i]
  geoexclude(xcomm=xcommused,yenv=yenvused,geodis=geodisused,distance=dis)
})
colnames(diss)<-c("10m","100m","500m","1000m","10000m","50000m","100000m","1000000m","5000000m")
save.file(diss,filename = "plant level_temp genus level_exclude",folder = save.wd)
stopCluster(c1)
