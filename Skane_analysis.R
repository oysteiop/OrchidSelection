######################################################################################
##### Analysis of trait-performance-fitness relationships in Orchids (2020-2021) #####
######################################################################################

# Analysis of Scanian populations

library(plyr)
library(lme4)
library(glmmTMB)
library(mgcv)
library(MuMIn)

rm(list=ls())

invlogit = function(x) 1/(1+exp(-x))

# Read data files ####
traitdata = read.table("data_2021/traitdata2021.txt", header=T)
fitnessdata = read.table("data_2021/fitnessdata2021.txt", header=T)

# Data formatting
traitdata$dataset = paste0(traitdata$site, "_", traitdata$year, "_", traitdata$period)
traitdata$dataset = as.factor(traitdata$dataset)
traitdata$ind = paste0(traitdata$dataset, "_", traitdata$individual)

fitnessdata$dataset = paste0(fitnessdata$site, "_", fitnessdata$year, "_", fitnessdata$period)
fitnessdata$dataset = as.factor(fitnessdata$dataset)
fitnessdata$ind = paste0(fitnessdata$dataset, "_", fitnessdata$individual)

# Insert missing species names
spVec = as.character(fitnessdata$species)
for(i in 1:nrow(fitnessdata)){
  if(is.na(fitnessdata$species[i])){
    spVec[i] = as.character(traitdata$species[which(traitdata$ind==fitnessdata$ind[i])[1]])
  }
}

fitnessdata$species = as.factor(spVec)

traitdata$dataset = paste0(traitdata$species, "_", traitdata$dataset)
traitdata$dataset = as.factor(traitdata$dataset)
traitdata$ind = paste0(traitdata$species, "_", traitdata$ind)

fitnessdata$dataset = paste0(fitnessdata$species, "_", fitnessdata$dataset)
fitnessdata$dataset = as.factor(fitnessdata$dataset)
fitnessdata$ind = paste0(fitnessdata$species, "_", fitnessdata$ind)

# Compile fitness data per individual ####
fitnessdata2 = subset(fitnessdata, select=c("ind", "flower", "pollinaria_removed", "pollen_on_stigma"))
fitnessdata2 = na.omit(fitnessdata2)

wdata = ddply(fitnessdata2, .(ind), summarize,
              n = sum(flower>0, na.rm=T),
              n_removed = sum(pollinaria_removed, na.rm=T),
              n_pollinated = sum(pollen_on_stigma, na.rm=T))

wdata$w_female = wdata$n_pollinated/wdata$n
wdata$w_male = wdata$n_removed/(wdata$n*2)
wdata$visited = 1*((wdata$n_removed+wdata$n_pollinated)>0) 

summary(wdata)
head(wdata)

# Mean traits per individual
traitmeans = ddply(traitdata, .(ind), summarize,
                   patch = unique(patch),
                   species = unique(species),
                   dataset = unique(dataset),
                   recorder_z = recorder[1],
                   height = mean(height_cm, na.rm=T),
                   flowers_open = mean(flowers_open, na.rm=T),
                   flowers_bud = mean(flowers_bud, na.rm=T),
                   flower_size = mean(corolla_height_mm, na.rm=T),
                   spur_length = mean(spur_length_mm, na.rm=T),
                   spur_width = mean(spur_width_mm, na.rm=T),
                   fruits_mean = mean(fruit_number, na.rm=T),
                   fruits_dev_mean = mean(fruits_small, na.rm=T))

traitmeans$fruits_total = apply(subset(traitmeans, select = c("fruits_mean", "fruits_dev_mean")), 1, sum, na.rm=T)
traitmeans$flowers = traitmeans$flowers_open + traitmeans$flowers_bud

head(traitmeans)

# Merge datafiles
dat = merge(traitmeans, wdata, by="ind", all=T)
dat = droplevels(dat)

names(dat)

fulldat = dat

#### Start of analyses ####

# Select dataset ####
levels(dat$dataset)

# All datasets
studies = levels(dat$dataset)[c(2,4,5,1,8)]
studies

betamat = matrix(NA, ncol=5, nrow=5)
SEmat = matrix(NA, ncol=5, nrow=5)
rownames(betamat) = rownames(SEmat) = studies

sslist = list()
prlist = list()
pmatlist = list()
rmflist = list()
mvlist = list()
mflist = list()
mmlist = list()
mwlist = list()
aiclist = list()
r2list = list()
mnaivelist = list()

nbot = 1000

for(s in 1:5){

seldata = studies[s]
seldata
dat = dat[dat$dataset==seldata, ]

# Removing outliers if any
if(seldata=="Dactylorhiza_majalis_Oerup_2021_NA"){
dat = dat[dat$height<50,]
dat = dat[dat$flower_size<16,]
}

# Centre and mean-scale
dat$height_c = scale(dat$height, scale=F)/mean(dat$height, na.rm=T)
dat$flowers_open_c = scale(dat$flowers_open, scale=F)/mean(dat$flowers_open, na.rm=T)
dat$flowers_c = scale(dat$flowers, scale=F)/mean(dat$flowers, na.rm=T)
dat$flower_size_c = scale(dat$flower_size, scale=F)/mean(dat$flower_size, na.rm=T)
dat$spur_length_c = scale(dat$spur_length, scale=F)/mean(dat$spur_length, na.rm=T)
dat$spur_width_c = scale(dat$spur_width, scale=F)/mean(dat$spur_width, na.rm=T)

# Summary statistics ####

# Trait means, variances, and coefficients of variation
summary_stats = data.frame(mean = apply(dat[,c("height", "flowers", "flower_size", "spur_length", "fruits_total", "visited", "w_female", "w_male")], 2, mean, na.rm=T),
                           sd = apply(dat[,c("height", "flowers", "flower_size", "spur_length", "fruits_total", "visited", "w_female", "w_male")], 2, sd, na.rm=T),
                           n_obs = apply(dat[,c("height", "flowers", "flower_size", "spur_length", "fruits_total",  "visited", "w_female", "w_male")], 2, function(x) sum(!is.na(x), na.rm=T)),
                           n_yes = apply(dat[,c("height", "flowers", "flower_size", "spur_length", "fruits_total",  "visited", "w_female", "w_male")], 2, function(x) sum(x>0, na.rm=T)))
summary_stats$cv = summary_stats$sd/summary_stats$mean
summary_stats$cv[6:8] = NA #The CV is not valid for the proportional variables

sslist[[s]] = summary_stats

# Phenotypic correlation matrix
P = cor(cbind(dat$height, dat$flowers, dat$flower_size, dat$spur_length), use="pairwise")
colnames(P) = rownames(P) = c("height", "flowers", "flower_size", "spur_length")

pmatlist[[s]] = P

# Fitting the component models ####

# Define visitation
rmf = cor(dat$w_female, dat$w_male, use="pairwise") # Correlation between pollen deposition and pollinarium removal

rmflist[[s]] = rmf

moddat = na.omit(subset(dat, select = c(patch, n, visited, w_male, w_female, 
                                        height_c, flowers, flowers_open_c, flower_size_c, spur_length_c)))

# Visitation
mv0 = glmmTMB(visited ~ 1 + (1|patch), family="binomial", data=moddat)

mv = glmmTMB(visited ~ height_c + flowers_open_c + flower_size_c + (1|patch), 
             family="binomial", data=moddat)


mvlist[[s]] = list(mv, signif(r.squaredGLMM(mv)[1,1], 3))

# Pollination conditional on visitation
vdat = moddat[which(moddat$visited>0),]

# Save pollination summary stats
prlist[[s]] = c(signif(sum(dat$visited, na.rm=T)/sum(dat$visited>-1, na.rm=T)*100, 3),
                round(mean(vdat$w_female, na.rm=T), 3),
                round(mean(vdat$w_male, na.rm=T), 3))

# Pollen deposition
mf0 = glmmTMB(w_female ~ 1 + (1|patch), family="binomial", weights=n, data=vdat)
mf = glmmTMB(w_female ~ height_c + flowers_open_c + flower_size_c + spur_length_c + (1|patch), family="binomial", weights=n, data=vdat)

mflist[[s]] = list(mf, signif(r.squaredGLMM(mf)[1,1], 3))

# Pollinarium removal
mm0 = glmmTMB(w_male ~ 1 + (1|patch), family="binomial", weights=n*2, data=vdat)
mm = glmmTMB(w_male ~ height_c + flowers_open_c + flower_size_c + spur_length_c + (1|patch), family="binomial", weights=n*2, data=vdat)

mmlist[[s]] = list(mm, signif(r.squaredGLMM(mm)[1,1], 3))

# Pollen to fruits analysis all data
wdatfull = na.omit(subset(fulldat, select = c(dataset, patch, n, visited, w_female, fruits_total, flowers)))
wdatfull = wdatfull[wdatfull$fruits_total>0,]

wdatfull$patch = as.factor(paste0(wdatfull$dataset, "_", wdatfull$patch))

mw0 = glmmTMB(fruits_total/flowers ~ 1 + (1|dataset/patch), "binomial", weights=flowers, dat=wdatfull)
mw = glmmTMB(fruits_total/flowers ~ w_female + (1|dataset/patch), "binomial", weights=flowers, dat=wdatfull)

if(sum(dat$fruits_total)>10){

# Population-specific fruit set analysis
wdat = na.omit(subset(dat, select = c(patch, n, visited, w_female, fruits_total, flowers)))
wdat = wdat[wdat$fruits_total>0,]

mw0 = glmmTMB(fruits_total/flowers ~ 1 + (1|patch), "binomial", weights=flowers, dat=wdat)
mw = glmmTMB(fruits_total/flowers ~ w_female + (1|patch), "binomial", weights=flowers, dat=wdat)

AIC(mw0, mw)

summary(mw)
}

mwlist[[s]] = list(mw, signif(r.squaredGLMM(mw)[1,1], 3))

# Compute selection gradients ####

# Visitation
rr = summary(mv)$coef$cond[,1]
pd = data.frame(Intercept=rep(1, nrow(dat)), subset(dat, select = names(rr[-1])))
predvis = invlogit(as.matrix(pd) %*% as.matrix(rr))

# Pollen deposition
rr = summary(mf)$coef$cond[,1]
pd = data.frame(Intercept=rep(1, nrow(dat)), subset(dat, select = names(rr[-1])))
predpoll = invlogit(as.matrix(pd) %*% as.matrix(rr))
predpoll = predpoll*predvis

# Fruitset
rr = summary(mw)$coef$cond[,1]
pd = data.frame(Intercept=rep(1, nrow(dat)), predpoll)
predfruitset = invlogit(as.matrix(pd) %*% as.matrix(rr))

# Fruits
predfruits = predfruitset*dat$flowers
dat$predfruits = predfruits

# Selection model
relfit = predfruits/mean(predfruits, na.rm=T)
ms = lm(relfit ~ flowers_open_c + height_c + flower_size_c + spur_length_c, dat=dat)
betas = summary(ms)$coef[,1]

betas*100

# Bootstrapping

betaList = list()

for(i in 1:nbot){
  
  # Visitation
  rr = MASS::mvrnorm(1, mu=summary(mv)$coef$cond[,1], Sigma=vcov(mv)[[1]])
  pd = data.frame(Intercept=rep(1, nrow(dat)), subset(dat, select = names(rr[-1])))
  predvis = invlogit(as.matrix(pd) %*% as.matrix(rr))
  
  # Pollen deposition
  rr = MASS::mvrnorm(1, mu=summary(mf)$coef$cond[,1], Sigma=vcov(mf)[[1]])
  pd = data.frame(Intercept=rep(1, nrow(dat)), subset(dat, select = names(rr[-1])))
  predpoll = invlogit(as.matrix(pd) %*% as.matrix(rr))
  
  predpoll = predpoll*predvis
  
  # Fruitset
  rr = MASS::mvrnorm(1, mu=summary(mw)$coef$cond[,1], Sigma=vcov(mw)[[1]])
  pd = data.frame(Intercept=rep(1, nrow(dat)), predpoll)
  predfruitset = invlogit(as.matrix(pd) %*% as.matrix(rr))
  
  # Fruits
  predfruits = predfruitset*dat$flowers
  
  # Selection model
  relfit = predfruits/mean(predfruits, na.rm=T)
  ms = lm(relfit~flowers_open_c + height_c + flower_size_c + spur_length_c, dat=dat)
  
  betaList[[i]] = summary(ms)$coef[,1]
  
}

betaList = lapply(betaList, function(x) as.data.frame(t(x)))
betaDat = rbind.fill(betaList)

apply(betaDat, 2, mean)*100
round(betas*100, 1)
round(apply(betaDat, 2, sd)*100, 2)
apply(betaDat, 2, quantile, c(0.025, 0.975))*100

betamat[s,] = betas*100
SEmat[s,] = apply(betaDat, 2, sd)*100

# Evaluating AIC support ####

AICff = AIC(mv) + AIC(mf)
AIC0 = AIC(mv0) + AIC(mf0)

aiclist[[s]] = AIC0-AICff

if(sum(dat$fruits_total)>10){
sel = which(dat$fruits_total>0)
r2 = cor(dat$fruits_total[sel], dat$predfruits[sel], "pairwise")^2 # Overall r^2 for fitness function

r2list[[s]] = r2

# Naive analysis
obs_relfit = dat$fruits_total/mean(dat$fruits_total, na.rm=T)
naivemod = lm(obs_relfit~height_c+flowers_open_c+flower_size_c+spur_length_c, na=na.exclude, dat=dat)

mnaivelist[[s]] = naivemod
}

dat = fulldat

}

# Save results
save(betamat, file="betamat.RData")
save(SEmat, file="SEmat.RData")

names(sslist) = names(prlist) = names(pmatlist) = names(rmflist) = names(mvlist) = 
names(mflist) = names(mmlist) = names(mwlist) = names(aiclist) = studies

names(r2list) = names(mnaivelist) = studies[1:3]

save(sslist, file="sslist.RData")
save(prlist, file="prlist.RData")
save(pmatlist, file="pmatlist.RData")
save(rmflist, file="rmflist.RData")
save(mvlist, file="mvlist.RData")
save(mflist, file="mflist.RData")
save(mmlist, file="mmlist.RData")
save(mwlist, file="mwlist.RData")
save(aiclist, file="aiclist.RData")
save(r2list, file="r2list.RData")
save(mnaivelist, file="mnaivelist.RData")

# Figures ####

# Visitation
coefs = summary(mv)$coef$cond
a1 = coefs[1,1]
b1_height = coefs[2,1]
b1_flowers = coefs[3,1]
b1_flower_size = coefs[4,1]

par(mfrow=c(1,3))

plot(dat$height_c, dat$visited)
xx = seq(min(dat$height_c, na.rm=T), max(dat$height_c, na.rm=T), 0.1)
yy = invlogit(a1 + b1_height*xx)
lines(xx,yy)

plot(dat$flowers_open_c, dat$visited)
xx = seq(min(dat$flowers_open_c, na.rm=T), max(dat$flowers_open_c, na.rm=T), 0.1)
yy = invlogit(a1 + b1_flowers*xx)
lines(xx,yy)

plot(dat$flower_size_c, dat$visited)
xx = seq(min(dat$flower_size_c, na.rm=T), max(dat$flower_size_c, na.rm=T), 0.1)
yy = invlogit(a1 + b1_flower_size*xx)
lines(xx,yy)

# Pollination
coefs = summary(mf)$coef$cond
a2 = coefs[1,1]
b2_height = coefs[2,1]
b2_flowers = coefs[3,1]
b2_flower_size = coefs[4,1]
b2_spur_length = coefs[5,1]

x11()
par(mfrow=c(2,2))

plot(vdat$height_c, vdat$w_female, las=1,
     xlab="Plant height",
     ylab="Proportion pollinated")
xx = seq(min(dat$height_c, na.rm=T), max(dat$height_c, na.rm=T), 0.01)
yy = invlogit(a2 + b2_height*xx)
lines(xx,yy)

plot(dat$flowers_open_c, dat$w_female, las=1,
     xlab="Open flowers",
     ylab="")
xx = seq(min(dat$flowers_open_c, na.rm=T), max(dat$flowers_open_c, na.rm=T), 0.01)
yy = invlogit(a2 + b2_flowers*xx)
lines(xx,yy)

plot(dat$flower_size_c, dat$w_female, las=1,
     xlab="Flower size",
     ylab="Proportion pollinated")
xx = seq(min(dat$flower_size_c, na.rm=T), max(dat$flower_size_c, na.rm=T), 0.01)
yy = invlogit(a2 + b2_flower_size*xx)
lines(xx,yy)

plot(dat$spur_length_c, dat$w_female, las=1,
     xlab="Spur length",
     ylab="")
xx = seq(min(dat$spur_length_c, na.rm=T), max(dat$spur_length_c, na.rm=T), 0.01)
yy = invlogit(a2 + b2_spur_length*xx)
lines(xx,yy)

# Fruitset
par(mfrow=c(1,1))
plot(wdat$w_female, wdat$fruits_total/wdat$flowers)

a3 = summary(mw)$coef[1,1]
b3_w_female = summary(mw)$coef[2,1]

xx = seq(min(dat$w_female, na.rm=T), max(dat$w_female, na.rm=T), 0.01)
yy = invlogit(a3 + b3_w_female*xx)
lines(xx, yy)

# Fitness
par(mfrow=c(1,1))

plot(dat$height, relfit)
plot(dat$flowers_open, relfit)
plot(dat$spur_length, relfit)
plot(dat$flower_size, relfit)


