


library(vegan)
library(parallel)
library(ecodist)
MRMnew<-function(formula = formula(data), data = sys.parent(), mrank = FALSE,standard.cc=FALSE,geoscale=c(0,Inf),nperm=999,no_cores=4,include.geo=TRUE){
  
  m <- match.call(expand.dots = FALSE)
  m2 <- match(c("formula", "data"), names(m), nomatch=0)
  m <- m[c(1, m2)]
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  m <- as.matrix(m)
  observed<-MRMonce(m=m,mrank =mrank, standard.cc=standard.cc,geoscale=geoscale,include.geo=include.geo)
  c1<-makeCluster(no_cores)
  permresult<-parSapply(c1,1:nperm,function(i){
    m2=m
    newY <- ecodist::full(m2[,1])
    newSample <- sample(nrow(newY))
    newY <- newY[newSample, newSample]
    m2[,1] <- ecodist::lower(newY)
    
    m2<-m2[m2[,2]>=geoscale[1] & m2[,2]< geoscale[2],]
    
    if(mrank) {
      m2 <- apply(m2, 2, rank)
    }
    
    
    if(standard.cc){
      m2<-scale(m2)
    }
    
    
    X <- m2[ ,2:ncol(m2), drop=FALSE]
    X <- cbind(rep(1, nrow(X)), X)
    Y <- m2[ ,1, drop=FALSE]
    
    if(include.geo){
      X<-X
    }else {X<-X[,-2]}
    
    
    
    nd <- nrow(X)
    
    # only need to calculate (X'X)^-1 once
    XX <- crossprod(X)
    XX <- solve(XX)
    
    # will need to calculate Xy for each permutation
    XY <- crossprod(X, Y)
    YY <- crossprod(Y)
    
    # regression coefficients
    b <- XX %*% XY
    rownames(b) <- c("Int", colnames(X)[2:ncol(X)])
    
    bXY <- crossprod(b, XY)
    SSE <- YY - bXY
    
    SSTO <- YY - sum(Y)^2/nd
    SSR = SSTO - SSE
    
    # R2 = 1 - SSE/SSTO
    R2 <- 1 - SSE/SSTO
    R2 <- as.vector(R2)
    
    # F* = MSR / MSE
    # MSR = SSR / (p - 1) 
    # MSE = SSE / (n - p)
    #b / sqrt(1 - R2)#
    p <- ncol(X) # number of parameters estimated
    F.value <- (SSR / (p - 1)) / (SSE / (nd - p))
    b.t<-b/sqrt(1-R2)
    result<-list(b.t=b.t,b=b,R2=R2,F.value=F.value)
    result
  },simplify = F)
 
   stopCluster(c1)
  ##combine permutations, generate a list of 4 elements##
  permresult2<-do.call(mapply, c(cbind, permresult))
  
  R2.all<-c(observed$R2,as.vector(t(permresult2$R2)))
  R2.pval <- length(R2.all[R2.all >= R2.all[1]])/(1+nperm)
  
  F.all <- c(observed$F.value,as.vector(t(permresult2$F.value)))
  F.pval <- length(F.all[F.all >= F.all[1]])/(1+nperm)
  
  # b.all contains pseudo-t of Legendre et al. 1994
  b.all<-t(cbind(observed$b.t,permresult2$b.t))
  b.pval <- apply(b.all, 2, function(x)length(x[abs(x) >= abs(x[1])])/(1+nperm))
  
  results <- list(coef=cbind(b=observed$b, pval=b.pval), r.squared=c(R2=observed$R2, pval = R2.pval),F.test=c(F.value=observed$F.value, F.pval = F.pval))
  results
  
}

MRMonce<-function(m, mrank = FALSE,standard.cc=FALSE,geoscale=c(0,Inf),include.geo=TRUE){
  
  m<-m[m[,2]>=geoscale[1] & m[,2]< geoscale[2],]
  
  if(mrank) {
    m <- apply(m, 2, rank)
  }
  
  
  if(standard.cc){
    m<-scale(m)
  }
  
  
  X <- m[ ,2:ncol(m), drop=FALSE]
  X <- cbind(rep(1, nrow(X)), X)
  Y <- m[ ,1, drop=FALSE]
  
  if(include.geo){
    X<-X
  }else {X<-X[,-2]}
  
  nd <- nrow(X)
  
  # only need to calculate (X'X)^-1 once
  XX <- crossprod(X)
  XX <- solve(XX)
  
  # will need to calculate Xy for each permutation
  XY <- crossprod(X, Y)
  YY <- crossprod(Y)
  
  # regression coefficients
  b <- XX %*% XY
  rownames(b) <- c("Int", colnames(X)[2:ncol(X)])
  
  bXY <- crossprod(b, XY)
  SSE <- YY - bXY
  
  SSTO <- YY - sum(Y)^2/nd
  SSR = SSTO - SSE
  
  # R2 = 1 - SSE/SSTO
  R2 <- 1 - SSE/SSTO
  R2 <- as.vector(R2)
  
  # F* = MSR / MSE
  # MSR = SSR / (p - 1) 
  # MSE = SSE / (n - p)
  #b / sqrt(1 - R2)#
  p <- ncol(X) # number of parameters estimated
  F.value <- (SSR / (p - 1)) / (SSE / (nd - p))
  b.t<-b/sqrt(1-R2)
  result<-list(b.t=b.t,b=b,R2=R2,F.value=F.value)
  result
}

##caution! the following code may contain many redundant ones..I haven't cleaned up yet

load("GW16S.taxadist.rda")
load("GW16S.unifrac.rda")


canberradis<-read.csv("canberra distance.csv",header = T,row.names = 1)
env.sample<-read.csv("resample.diversity and env_sample level.csv",row.names = 1,header = T)

env.sample<-env.sample[match(row.names(canberradis),row.names(env.sample)),]

env<-data.frame(temp=env.sample$Air.temperature.Mean.annual,precip=env.sample$Precipitation.Annual,BOD=env.sample$Influent.BOD..mg.L.,
                Sludge.age=env.sample$Sludge.Age..days.,pH=env.sample$pH,DO=env.sample$DO)

row.names(env)<-row.names(env.sample)
geo<-env.sample[,17:18]
source("distLLE.r")
geodis=distLLE(latitude = geo$Latitude.Google,longitude = geo$Longitude.Google,site.name = rownames(geo),short = FALSE)
min.dist=2 #m
lngeodist=log(geodis+min.dist)


#global
library(vegan)
idused<-complete.cases(env)
lngeodisused<-lngeodist[idused,idused]
commdisused<-commdist[idused,idused]
envused<-env[idused,]
tempdis<-vegdist(scale(envused$temp),method="euclidean")
precipdis<-vegdist(scale(envused$precip),method="euclidean")
BODdis<-vegdist(scale(envused$BOD),method="euclidean")
Sludgedis<-vegdist(scale(envused$Sludge.age),method="euclidean")
pHdis<-vegdist(scale(envused$pH),method="euclidean")
DOdis<-vegdist(scale(envused$DO),method="euclidean")


library(ecodist)

MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused)+tempdis+precipdis+BODdis+Sludgedis+pHdis+DOdis,standard.cc = T,no_cores=5)


MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused)+tempdis+BODdis+Sludgedis+pHdis+DOdis,standard.cc = T,no_cores=5)

MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused),standard.cc = T)
MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused)+tempdis+BODdis+Sludgedis+pHdis+DOdis,standard.cc = T,no_cores=5,include.geo = F)


##within city
Commdist<-as.dist(log(commdisused+0.005))
geodist<-as.dist(lngeodisused)
scale=c(0,log(1e5+2))
MRMnew(Commdist~geodist+tempdis+precipdis+BODdis+Sludgedis+pHdis+DOdis,standard.cc = T,geoscale =scale)

env<-data.frame(
                Sludge.age=env.sample$Sludge.Age..days.)

idused<-complete.cases(env)
lngeodisused<-lngeodist[idused,idused]
commdisused<-commdist[idused,idused]
envused<-env[idused,]


Sludgedis<-vegdist(scale(envused),method="euclidean")

MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused)+Sludgedis,standard.cc = T,no_cores=5,geoscale =scale,include.geo = T)
MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused),standard.cc = T,no_cores=5,geoscale =scale,include.geo = T)
MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused)+Sludgedis,standard.cc = T,no_cores=5,geoscale =scale,include.geo = F)

##across city, within a continent
env<-data.frame(temp=env.sample$Air.temperature.Mean.annual,precip=env.sample$Precipitation.Annual,BOD=env.sample$Influent.BOD..mg.L.,
                Sludge.age=env.sample$Sludge.Age..days.,pH=env.sample$pH,DO=env.sample$DO)

idused<-complete.cases(env)
lngeodisused<-lngeodist[idused,idused]
commdisused<-commdist[idused,idused]
envused<-env[idused,]
tempdis<-vegdist(scale(envused$temp),method="euclidean")
precipdis<-vegdist(scale(envused$precip),method="euclidean")
BODdis<-vegdist(scale(envused$BOD),method="euclidean")
Sludgedis<-vegdist(scale(envused$Sludge.age),method="euclidean")
pHdis<-vegdist(scale(envused$pH),method="euclidean")
DOdis<-vegdist(scale(envused$DO),method="euclidean")

scale=c(log(1e5+2),log(5e6+2))
MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused)+tempdis+precipdis+BODdis+Sludgedis+pHdis+DOdis,standard.cc = T,no_cores=5,geoscale =scale,include.geo = T)


env<-data.frame(temp=env.sample$Air.temperature.Mean.annual,
                Sludge.age=env.sample$Sludge.Age..days.,DO<-env.sample$DO,pH=env.sample$pH)



idused<-complete.cases(env)
lngeodisused<-lngeodist[idused,idused]
commdisused<-commdist[idused,idused]
envused<-env[idused,]

tempdis<-vegdist(scale(envused$temp),method="euclidean")
#precipdis<-vegdist(scale(envused$precip),method="euclidean")
Sludgedis<-vegdist(scale(envused$Sludge.age),method="euclidean")
pHdis<-vegdist(scale(envused$pH),method="euclidean")
DOdis<-vegdist(scale(envused$DO),method="euclidean")
MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused)+tempdis+Sludgedis+DOdis+pHdis,standard.cc = T,no_cores=5,geoscale =scale,include.geo = T)


MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused)+tempdis+precipdis+Sludgedis+pHdis,standard.cc = T,no_cores=5,geoscale =scale,include.geo = T)
MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused)+tempdis+precipdis+Sludgedis+pHdis,standard.cc = T,no_cores=5,geoscale =scale,include.geo = F)
MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused),standard.cc = T,no_cores=5,geoscale =scale,include.geo = T)


##across continent##
env<-data.frame(temp=env.sample$Air.temperature.Mean.annual,BOD=env.sample$Influent.BOD..mg.L.,
                Sludge.age=env.sample$Sludge.Age..days.,pH=env.sample$pH,DO=env.sample$DO)

idused<-complete.cases(env)
lngeodisused<-lngeodist[idused,idused]
commdisused<-commdist[idused,idused]
envused<-env[idused,]
tempdis<-vegdist(scale(envused$temp),method="euclidean")

BODdis<-vegdist(scale(envused$BOD),method="euclidean")
Sludgedis<-vegdist(scale(envused$Sludge.age),method="euclidean")
pHdis<-vegdist(scale(envused$pH),method="euclidean")
DOdis<-vegdist(scale(envused$DO),method="euclidean")

scale=c(log(5e6+2),Inf)
MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused)+tempdis+precipdis+BODdis+Sludgedis+pHdis+DOdis,standard.cc = T,no_cores=5,geoscale =scale,include.geo = T)
MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused)+tempdis+BODdis+Sludgedis+pHdis+DOdis,standard.cc = T,no_cores=5,geoscale =scale,include.geo = T)


MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused)+tempdis+BODdis+Sludgedis+pHdis+DOdis,standard.cc = T,no_cores=5,geoscale =scale,include.geo = F)
MRMnew(as.dist(log(commdisused+0.005))~as.dist(lngeodisused),standard.cc = T,no_cores=5,geoscale =scale,include.geo = T)


