
##richness and latitude##
dat<-read.table("latitude diversity for plot.csv",sep=",",row.names=1,header=T)
hemi<-ifelse(dat$Latitude.Google>0,"North","South")
dat2<-data.frame(cbind(dat,hemi=hemi))
library(ggplot2)
dat2$Temperature
ggplot(dat2,aes(x=lat.abs,y=phylo0))+
  geom_point(size=2,alpha=0.6,aes(shape=hemi,colour=Temperature))+
  geom_smooth(method = "lm",formula =y ~ poly(x, 2))+
  # scale_colour_manual(values = c("North"="#ff0000","South"="#0000ff"))+
  scale_shape_manual(values = c(19,0))+
  scale_colour_gradientn(colours =c("#aa0000ff","#d40000ff","#ff0000ff","#ff5555ff","#ff8080ff","#ffaaaaff"),trans = 'reverse')+
  xlab("Absolute latitude") +
  ylab("Phylogenetic diversity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

summary(lm(phylo0~I(lat.abs)+I(lat.abs^2), data=dat2))

##average copy number and BOD/1+recycling ratio##
dat4=data.frame(copy=dat$average.copy.number,BODratio=dat$BOD.recycling.ratio,BOD=dat$Influent.BOD..mg.L.)

ggplot(dat4,aes(x=BODratio,y=copy))+
  geom_point(size=3,alpha=0.5,colour="#0066ffff")+
  geom_smooth(method = "lm",colour="red")+
  xlab("BOD/(1+recycling ratio)") +
  ylab("Average rRNA gene copy number") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

summary(lm(copy~BODratio, data=dat4))

