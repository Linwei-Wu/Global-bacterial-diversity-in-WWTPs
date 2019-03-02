

library(lavaan)
dat<-read.csv("ENV and comm data For SEM.csv",row.names = 1,header=T)

dat$tempSQ=dat$Air.temperature.Mean.annual^2
dat$SludgeageSQ=dat$Sludge.Age^2

##scale to the same variance
dat<-data.frame(scale(dat,center = F))
pairs(dat)



##model
model <- '
# composite variable

Sludge.effect<~ 0.4*Sludge.Age+SludgeageSQ

# regressions

MLSS~F.M+Air.temperature.Mean.annual+DO+Sludge.Age

Sludge.Age~F.M+Influent.BOD+Air.temperature.Mean.annual
F.M~Sludge.Age+Influent.BOD
SludgeageSQ~Influent.BOD
S~Sludge.effect+Influent.BOD+MLSS+DO+PC1.bray+Air.temperature.Mean.annual
PC1.bray~Sludge.effect+Air.temperature.Mean.annual+F.M+Influent.BOD+MLSS+DO
BOD.removal~MLSS+PC1.bray+F.M+DO+Sludge.effect
#NH4.removal~Influent.BOD+PC1.phy+F.M
#TN.removal~Influent.BOD+PC1.phy+F.M+DO

# residual correlations
SludgeageSQ ~~ Sludge.Age

'
abioticCompFit <- sem(model, missing="ml",data=dat)
summary(abioticCompFit, rsquare=T, standardized=T)
summary(abioticCompFit, fit.measures=TRUE)
residuals(abioticCompFit, type="cor")
modificationIndices(abioticCompFit,standardized=F)

