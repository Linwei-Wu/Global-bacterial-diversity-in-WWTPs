

otu.file="~/16S.OTU.before.Resample/GwMC.16SUP.otu.uc.all.txt"
prefix="GW16S.uc.rmsg.plant"
#prefix="GW16S.uc.rmsg.site"
out.wd="~/rarefaction"
sampcate.file="~Sample.Continent.Country.Site.WWTP.csv"



#######################################
# internal use. not published yet. please at least acknowledge if you used.
#install.packages("ieggr",repos = c("ftp://gr:gr!123@129.15.40.254","http://cran.us.r-project.org"))
# In some cases, you need to install it in R rather than RStudio.

library(ieggr)
library(vegan)
out.wd=sub("file:///","",out.wd)
setwd(out.wd)

com=lazyopen(otu.file)
dim(com)
comm=t(com)


samp.grp=lazyopen(sampcate.file)

##remove global singleton
comm<-comm[match(row.names(samp.grp),row.names(comm)),]
sum(rowSums(comm)==0);sum(colSums(comm)==0)


Comm<-comm[,colSums(comm)>1]
sum(colSums(Comm)<=1)

samp.grp=data.frame(samples=row.names(samp.grp),samp.grp)
samp.group=samp.grp[,8,drop=FALSE]

group.rcv<-function(comm,samp.group,prefix=character(0),geometric=TRUE,ser.num=500,out.wd=getwd(),save.out=TRUE)
{
  library(ieggr)
  out.wd=sub("file:///","",out.wd)
  grp.lev=levels(as.factor(samp.group[,1]))
  grp.num=length(grp.lev)
  
  coma=sapply(1:grp.num,
              function(i)
              {
                samps=rownames(samp.group)[which(samp.group[,1]==grp.lev[i])]
                colSums(comm[match(samps,rownames(comm)),])
              })
  #dim(coma)
  colnames(coma)=grp.lev
  
  serf<-function(ser.num,J,geometric=geometric)
  {
    if(geometric)
    {
      ser=unique(round(exp((log(J)/ser.num)*(1:ser.num))))
      if(ser[1]!=1){ser=c(1,ser)}
      if(ser[length(ser)]!=J){ser=c(ser,J)}
    }else{
      ser=floor(J/ser.num)*(1:ser.num)
      if(ser[1]!=1){ser=c(1,ser)}
      if(ser[length(ser)]!=J){ser=c(ser,(J-1),J)}
    }
    ser=c(ser[1:(length(ser)-1)],ser[length(ser)]-1,ser[length(ser)])
    ser
  }
  
  rare=list()
  rslop=list()
  
  for(i in 1:grp.num)
  {
    message("group i=",i," in ",grp.num,". ",date())
    seri=serf(ser.num = ser.num,J=sum(coma[,i]),geometric = geometric)
    rarei=rarefy(coma[,i,drop=FALSE],sample = seri,MARGIN = 2)
    rare[[i]]=data.frame(sequence.number=as.vector(seri),OTU.number=as.vector(rarei))
    rslopi=rareslope(t(coma[,i,drop=FALSE]),sample = seri)
    rslop[[i]]=data.frame(sequence.number=as.vector(seri),rare.slope=as.vector(rslopi))
  }
  names(rare)<-names(rslop)<-grp.lev
  if(save.out)
  {
    save.file(coma,prefix = prefix,filename = "group.AbSum",folder = out.wd)
    
    list.cbind<-function(lista)
    {
      len=length(lista)
      rownum=max(sapply(1:len,function(i){nrow(lista[[i]])}))
      if(is.null(names(lista))){names(lista)=paste("l",1:len,sep = "")}
      out.list=lapply(1:len,
                      function(i)
                      {
                        comple=matrix("",nrow = (rownum-nrow(lista[[i]])),ncol=ncol(lista[[i]]))
                        outl=rbind(as.matrix(lista[[i]]),comple)
                        if(is.null(colnames(lista[[i]]))){colnames(lista[[i]])=1:ncol(lista[[i]])}
                        colnames(outl)=paste(names(lista)[i],colnames(lista[[i]]),sep = ".")
                        outl
                      })
      Reduce(cbind,out.list)
    }
    rare.all=list.cbind(rare)
    rslop.all=list.cbind(rslop)
    save.file(rare.all,prefix = prefix,filename = "raref.curve",folder = out.wd)
    save.file(rslop.all,prefix = prefix,filename = "raref.slope",folder = out.wd)
  }
  list(group.sum=coma,rarefy=rare,rare.slope=rslop)
}

######sample rarefaction curve
sample.rcv<-function(comm,prefix=character(0),geometric=TRUE,ser.num=500,out.wd=getwd(),save.out=TRUE)
{
  library(ieggr)
  out.wd=sub("file:///","",out.wd)
  grp.lev=row.names(comm)
  grp.num=length(grp.lev)
  
  coma=t(comm)
  #dim(coma)
  colnames(coma)=grp.lev
  
  serf<-function(ser.num,J,geometric=geometric)
  {
    if(geometric)
    {
      ser=unique(round(exp((log(J)/ser.num)*(1:ser.num))))
      if(ser[1]!=1){ser=c(1,ser)}
      if(ser[length(ser)]!=J){ser=c(ser,J)}
    }else{
      ser=floor(J/ser.num)*(1:ser.num)
      if(ser[1]!=1){ser=c(1,ser)}
      if(ser[length(ser)]!=J){ser=c(ser,(J-1),J)}
    }
    ser=c(ser[1:(length(ser)-1)],ser[length(ser)]-1,ser[length(ser)])
    ser
  }
  
  rare=list()
  rslop=list()
  
  for(i in 1:grp.num)
  {
    message("group i=",i," in ",grp.num,". ",date())
    seri=serf(ser.num = ser.num,J=sum(coma[,i]),geometric = geometric)
    rarei=rarefy(coma[,i,drop=FALSE],sample = seri,MARGIN = 2)
    rare[[i]]=data.frame(sequence.number=as.vector(seri),OTU.number=as.vector(rarei))
    rslopi=rareslope(t(coma[,i,drop=FALSE]),sample = seri)
    rslop[[i]]=data.frame(sequence.number=as.vector(seri),rare.slope=as.vector(rslopi))
  }
  names(rare)<-names(rslop)<-grp.lev
  if(save.out)
  {
    save.file(coma,prefix = prefix,filename = "group.AbSum",folder = out.wd)
    
    list.cbind<-function(lista)
    {
      len=length(lista)
      rownum=max(sapply(1:len,function(i){nrow(lista[[i]])}))
      if(is.null(names(lista))){names(lista)=paste("l",1:len,sep = "")}
      out.list=lapply(1:len,
                      function(i)
                      {
                        comple=matrix("",nrow = (rownum-nrow(lista[[i]])),ncol=ncol(lista[[i]]))
                        outl=rbind(as.matrix(lista[[i]]),comple)
                        if(is.null(colnames(lista[[i]]))){colnames(lista[[i]])=1:ncol(lista[[i]])}
                        colnames(outl)=paste(names(lista)[i],colnames(lista[[i]]),sep = ".")
                        outl
                      })
      Reduce(cbind,out.list)
    }
    rare.all=list.cbind(rare)
    rslop.all=list.cbind(rslop)
    save.file(rare.all,prefix = prefix,filename = "raref.curve",folder = out.wd)
    save.file(rslop.all,prefix = prefix,filename = "raref.slope",folder = out.wd)
  }
  list(group.sum=coma,rarefy=rare,rare.slope=rslop)
}

#######################################################################
grp.rcv=group.rcv(comm = Comm,samp.group=samp.group,prefix=prefix,geometric=TRUE,ser.num=500,out.wd=out.wd,save.out=TRUE)

sample.rcv2=sample.rcv(comm = Comm,prefix=prefix,geometric=TRUE,ser.num=500,out.wd=out.wd,save.out=TRUE)



curvefile<-read.csv("GW16S.uc.rmsg.plant.raref.curve.csv",header = T,row.names = 1)
names(curvefile)

#x1<-sub("(\\..*?)\\.","\\1",x)
#sub("(.*?)(\\.)(.*)","\\3\\2\\1",x1)

name1<-sub("(\\..*?)\\.","\\1",names(curvefile)[1:8])
name2<-sub("(\\..*?)\\.","\\1",sub("\\.","",names(curvefile)[9:12]))
name<-sub("(.*?)(\\.)(.*)","\\3\\2\\1",c(name1,name2))

name<-sub("(.*?)(\\.)(.*)","\\3\\2\\1",name1)

names(curvefile)<-name
curve2<-reshape(curvefile, direction="long", varying=TRUE, sep=".", timevar="Plant")
curvedata<-curve2[!is.na(curve2$sequencenumber),]
head(curvedata)
conti<-sapply(1:nrow(curvedata),function(i){
  samp.grp[curvedata$Sample[i],2]})

conti<-samp.grp2$Continent[match(curvedata$Plant,samp.grp2$WWTPID)]

curvedata=data.frame(curvedata,Continent=conti)
library(RColorBrewer)
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n = 6, name = "Dark2")

library(ggplot2)
ggplot(data = curvedata, aes(x = sequencenumber, y = OTUnumber,group=Sample,  color = Continent))+
  geom_line()+
  scale_colour_manual(values=brewer.pal(n = 6, name = "Dark2"))+
  labs(x = "Sequence number", y = "OTU number")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(data = curvedata, aes(x = sequencenumber, y = OTUnumber,group=Plant,  color = Continent))+
  geom_line()+
  scale_colour_manual(values=brewer.pal(n = 6, name = "Dark2"))+
  labs(x = "Sequence number", y = "OTU number")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

