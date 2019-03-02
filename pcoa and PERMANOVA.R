
#read otu table#

com.file="otu table170612.csv"
tree.file="GW16SUP.ucrmsg.res.tree.nwk"


# internal use. not published yet. please at least acknowledge if you used.
#install.packages("ieggr",repos = c("ftp://gr:gr!123@129.15.40.254","http://cran.us.r-project.org"))
# In some cases, you need to install it in R rather than RStudio.
library(ieggr)
library(ape)
wd=iwd(wd)

#taxonomic 
comm=t(lazyopen(com.file))
beta.binary=beta.g(log(comm+1),dist.method = c("jaccard","bray"),abundance.weighted = FALSE)
beta.weight=beta.g(log(comm+1),dist.method = c("jaccard","bray"),abundance.weighted = TRUE)

betaTD=list(jaccard=beta.binary$jaccard.uw,ruzicka=beta.weight$ruzicka.wt,sorensen=beta.binary$sorensen.uw,bray=beta.weight$bray.wt)
save(betaTD,file="GW16S.taxadist.rda")


#phylogenetic
tree=lazyopen(tree.file)

source("unifrac.g.r")
tt=Sys.time()
ufn=unifrac.g(otu.tab =log(comm+1),tree = tree,nworker=4)
(t1=format(Sys.time()-tt))
save(ufn,file="GW16S.unifrac.rda")


#######################################################
grpin=read.csv("Sample.categorical env.csv",row.names = 1,header=T)

load("GW16S.taxadist.rda")

betau<-as.matrix(betaTD$bray)

pco=pcoa(betau)
eigen=pco$values
samp.t3=pco$vectors[,1:3]

grpin2<-grpin[match(row.names(samp.t3),row.names(grpin)),]
sum(row.names(samp.t3)==row.names(grpin2))

data<-data.frame(cbind(samp.t3,grpin2))

##ggplot2##
library(ggplot2)
data$Continent
ggplot(data,aes(x=Axis.1,y=Axis.2,colour=Continent))+
  geom_point()+
  xlab("pcoa 1")+
  ylab("pcoa 2")


data2=data[!is.na(data$Climate.type),]
ggplot(data2,aes(x=Axis.1,y=Axis.2,colour=Climate.type))+
  geom_point()+
  xlab("pcoa 1")+
  ylab("pcoa 2")





##adonis test##
load("GW16S.unifrac.rda")
load("GW16S.taxadist.rda")

grp<-read.csv("metadata for R.csv", header=T, row.names = 1)

braydist<-as.matrix(betaTD$bray)
braydist<-braydist[match(row.names(grp),row.names(braydist)),match(row.names(grp),colnames(braydist))]
sum(row.names(grp)==row.names(braydist))


unif<-ufn$d_1
unif<-unif[match(row.names(grp),row.names(unif)),match(row.names(grp),colnames(unif))]
sum(row.names(grp)==row.names(unif))


##Africa vs Asia
levels(as.factor(grp2$Continent))

id<-which(grp2$Continent %in% c("Africa","Asia"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])

##Africa vs Australasia
id<-which(grp2$Continent %in% c("Africa","Australasia"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])


##Africa vs Europe
id<-which(grp2$Continent %in% c("Africa","Europe"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])

##Africa vs North America
id<-which(grp2$Continent %in% c("Africa","North America"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])

##Africa vs South America
id<-which(grp2$Continent %in% c("Africa","South America"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])

##Asia vs Australasia
id<-which(grp2$Continent %in% c("Asia","Australasia"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])

##Asia vs Europe
id<-which(grp2$Continent %in% c("Asia","Europe"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])

##Asia vs North America
id<-which(grp2$Continent %in% c("Asia","North America"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])

##Asia vs South America
id<-which(grp2$Continent %in% c("Asia","South America"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])


##Australasia vs Europe
id<-which(grp2$Continent %in% c("Australasia","Europe"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])


##Australasia vs North America
id<-which(grp2$Continent %in% c("Australasia","North America"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])

##Australasia vs South America
id<-which(grp2$Continent %in% c("Australasia","South America"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])

##Europe vs North America
id<-which(grp2$Continent %in% c("Europe","North America"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])

##Europe vs North America
id<-which(grp2$Continent %in% c("Europe","North America"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])

##Europe vs South America
id<-which(grp2$Continent %in% c("Europe","South America"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])

##North America vs South America
id<-which(grp2$Continent %in% c("North America","South America"))
adonis(as.dist(braydist[id,id])~Continent, data=grp2[id,])


###adonis of continent, climate type and activated sludge

braydist<-betaTD$bray

unif<-ufn$d_1
sum(row.names(grp2)==row.names(unif))
sum(row.names(grp2)==row.names(braydist))

grp2<-grp
grp2$Climate.type2<-substring(grp2$Climate.type, 1, 1)
grp2$Climate.type2<-as.factor(grp2$Climate.type2)

id<-which(!is.na(grp2$activated.sludge.type))

adonis(as.dist(braydist[id,id])~Continent*Climate.type2*activated.sludge.type,data=grp2[id,])


