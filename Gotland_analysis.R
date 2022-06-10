######################################################################################
##### Analysis of trait-performance-fitness relationships in Orchids (2020-2021) #####
######################################################################################

# Analysis of Gotland population

rm(list=ls())

library(plyr)
library(lme4)
library(glmmTMB)
library(mgcv)
library(MuMIn)

invlogit = function(x) 1/(1+exp(-x))

# Read data
datadir = "data_CRO"

traitdata <- read.csv(paste0(datadir, "/Gotland_trait_data.csv"), sep=";", dec=",")
fitnessdata <- read.csv(paste0(datadir, "/Gotland_fitness_data.csv"),  sep=";", dec=",")

# Compile fitness data per individual
fitnessdata = na.omit(fitnessdata)

wdata = ddply(fitnessdata, .(individual), summarise,
              n = sum(flower>0, na.rm=T),
              n_removed = sum(2-pollinaria_remaining, na.rm = T),
              n_pollinated = sum(pollen_on_stigma, na.rm = T))
wdata$w_female = wdata$n_pollinated/wdata$n
wdata$w_male = wdata$n_removed/(wdata$n*2) ## Two pollinaria per flower

summary(wdata)
head(wdata)

rmf = cor(wdata$w_female, wdata$w_male)
Gotland_rmf = save(rmf, file="Gotland_rmf.RData")

# Define visitation
visited <- 1*((wdata$w_female + wdata$w_male)>0)
wdata$visited <- visited
sum(visited, na.rm = TRUE)

# Mean traits per individual
names(traitdata)

traitmeans = ddply(traitdata, .(Site, Period, Species, individual), summarise,
                   patch = unique(Patch),
                   height = mean(height_cm, na.rm=T),
                   flowers_open = mean(flowers_open, na.rm=T),
                   flowers_bud = mean(flowers_bud, na.rm=T),
                   flower_size = mean(flower_size_mm, na.rm=T),
                   spur_length = mean(spur_length_mm, na.rm=T),
                   spur_width = mean(spur_width_mm, na.rm=T),
                   fruit_number= mean(fruit_number, na.rm=T))
head(traitmeans)

# Merge datafiles
dat = merge(wdata, traitmeans, by="individual", all=T)
names(dat)

dat$flowers = dat$flowers_open + dat$flowers_bud

# Check for outliers
hist(dat$height)
hist(dat$flowers_open)
hist(dat$flower_size)
hist(dat$spur_length)
hist(dat$spur_width)

# Centre and mean-scale
dat$height_c = scale(dat$height, scale=F)/mean(dat$height, na.rm=T)
dat$flowers_open_c = scale(dat$flowers_open, scale=F)/mean(dat$flowers_open, na.rm=T)
dat$flowers_c = scale(dat$flowers, scale=F)/mean(dat$flowers, na.rm=T)
dat$flower_size_c = scale(dat$flower_size, scale=F)/mean(dat$flower_size, na.rm=T)
dat$spur_length_c = scale(dat$spur_length, scale=F)/mean(dat$spur_length, na.rm=T)
dat$spur_width_c = scale(dat$spur_width, scale=F)/mean(dat$spur_width, na.rm=T)

# Summary statistics ####

# Trait means, variances, and coefficients of variation
summary_stats = data.frame(mean = apply(dat[,c("height", "flowers", "flower_size", "spur_length", "fruit_number", "visited", "w_female", "w_male")], 2, mean, na.rm=T),
                           sd = apply(dat[,c("height", "flowers", "flower_size", "spur_length", "fruit_number", "visited", "w_female", "w_male")], 2, sd, na.rm=T),
                           n_obs = apply(dat[,c("height", "flowers", "flower_size", "spur_length", "fruit_number",  "visited", "w_female", "w_male")], 2, function(x) sum(!is.na(x), na.rm=T)),
                           n_yes = apply(dat[,c("height", "flowers", "flower_size", "spur_length", "fruit_number",  "visited", "w_female", "w_male")], 2, function(x) sum(x>0, na.rm=T)))
summary_stats$cv = summary_stats$sd/summary_stats$mean
summary_stats$cv[6:8] = NA #The CV is not valid for the proportional variables

signif(summary_stats, 4)

Gotland_ss = save(summary_stats, file="Gotland_ss.RData")

# Phenotypic correlation matrix ####
P = cor(cbind(dat$height, dat$flowers, dat$flower_size, dat$spur_length,dat$spur_width), use="pairwise")
colnames(P) = rownames(P) = c("height", "flowers", "flower_size", "spur_length", "spur_width")
signif(P, 3)

Gotland_pmat = save(P, file="Gotland_pmat.RData")

# Fitting the component models ####

# Define visitation
dat$visited = 1*((dat$n_removed+dat$n_pollinated)>0) 

sum(dat$visited, na.rm=T)
sum(dat$visited, na.rm=T)/sum(dat$visited>-1, na.rm=T)

moddat = na.omit(subset(dat, select = c(n, visited, w_male, w_female, 
                                        height_c, flowers, flowers_open_c, flower_size_c, spur_length_c, spur_width_c)))

mv0 = glm(visited ~ 1, family="binomial", data=moddat)

mv = glm(visited ~ height_c + flowers_open_c + flower_size_c, 
             family="binomial", data=moddat)

mv1 = glm(visited ~ height_c + flowers_open_c, 
              family="binomial", data=moddat)
mv2 = glm(visited ~ height_c + flower_size_c, 
              family="binomial", data=moddat)
mv3 = glm(visited ~ flowers_open_c + flower_size_c, 
              family="binomial", data=moddat)
mv4 = glm(visited ~ height_c, 
              family="binomial", data=moddat)
mv5 = glm(visited ~ flowers_open_c, 
              family="binomial", data=moddat)
mv6 = glm(visited ~ flower_size_c, 
              family="binomial", data=moddat)

AICtab = AIC(mv0, mv, mv1, mv2, mv3, mv4, mv5, mv6)
AICtab[order(AICtab$AIC),]

mvsum = list(mv, signif(r.squaredGLMM(mv)[1,1], 3))

Gotland_mv = save(mvsum, file="Gotland_mv.RData")

# Pollination conditional on visitation
vdat = moddat[which(moddat$visited>0),]

percentage_visited = (sum(dat$visited, na.rm = T)/sum(!is.na(dat$visited)))*100
percentage_visited
mean(vdat$w_female, na.rm=T)
mean(vdat$w_male, na.rm=T)

Gotland_pr = c(percentage_visited, mean(vdat$w_female), mean(vdat$w_male))
save(Gotland_pr, file = "Gotland_pr.RData")

# Pollen deposition
mf0 = glm(w_female ~ 1, family="binomial", weights=n, data=vdat)
mf = glm(w_female ~ height_c + flowers_open_c + flower_size_c + spur_length_c, family="binomial", weights=n, data=vdat)

mf1 = glm(w_female ~ height_c + flowers_open_c + flower_size_c, family="binomial", weights=n, data=vdat)
mf2 = glm(w_female ~ height_c + flowers_open_c + spur_length_c, family="binomial", weights=n, data=vdat)
mf3 = glm(w_female ~ height_c + flower_size_c + spur_length_c, family="binomial", weights=n, data=vdat)
mf4 = glm(w_female ~ flowers_open_c + flower_size_c + spur_length_c, family="binomial", weights=n, data=vdat)
mf5 = glm(w_female ~ height_c + flowers_open_c, family="binomial", weights=n, data=vdat)
mf6 = glm(w_female ~ height_c + flower_size_c, family="binomial", weights=n, data=vdat)
mf7 = glm(w_female ~ height_c + spur_length_c, family="binomial", weights=n, data=vdat)
mf8 = glm(w_female ~ flowers_open_c + flower_size_c, family="binomial", weights=n, data=vdat)
mf9 = glm(w_female ~ flowers_open_c + spur_length_c, family="binomial", weights=n, data=vdat)
mf10 = glm(w_female ~ flower_size_c + spur_length_c, family="binomial", weights=n, data=vdat)
mf11 = glm(w_female ~ height_c, family="binomial", weights=n, data=vdat)
mf12 = glm(w_female ~ flowers_open_c, family="binomial", weights=n, data=vdat)
mf13 = glm(w_female ~ flower_size_c, family="binomial", weights=n, data=vdat)
mf14 = glm(w_female ~ spur_length_c, family="binomial", weights=n, data=vdat)

AICtab = AIC(mf0, mf, mf1, mf2, mf3, mf4, mf5, mf6, mf7, mf8, mf9, mf10, mf11, mf12, mf13, mf14)
AICtab[order(AICtab$AIC),]

mfsum = list(mf, signif(r.squaredGLMM(mf)[1,1], 3))
Gotland_mf = save(mfsum, file="Gotland_mf.RData")

# Pollinarium removal
mm0 = glm(w_male ~ 1, family="binomial", weights=n*2, data=vdat)
mm = glm(w_male ~ height_c + flowers_open_c + flower_size_c + spur_length_c, family="binomial", weights=n*2, data=vdat)

mm1 = glm(w_male ~ height_c + flowers_open_c + flower_size_c, family="binomial", weights=n*2, data=vdat)
mm2 = glm(w_male ~ height_c + flowers_open_c + spur_length_c, family="binomial", weights=n*2, data=vdat)
mm3 = glm(w_male ~ height_c + flower_size_c + spur_length_c, family="binomial", weights=n*2, data=vdat)
mm4 = glm(w_male ~ flowers_open_c + flower_size_c + spur_length_c, family="binomial", weights=n*2, data=vdat)
mm5 = glm(w_male ~ height_c + flowers_open_c, family="binomial", weights=n*2, data=vdat)
mm6 = glm(w_male ~ height_c + flower_size_c, family="binomial", weights=n*2, data=vdat)
mm7 = glm(w_male ~ height_c + spur_length_c, family="binomial", weights=n*2, data=vdat)
mm8 = glm(w_male ~ flowers_open_c + flower_size_c, family="binomial", weights=n*2, data=vdat)
mm9 = glm(w_male ~ flowers_open_c + spur_length_c, family="binomial", weights=n*2, data=vdat)
mm10 = glm(w_male ~ flower_size_c + spur_length_c, family="binomial", weights=n*2, data=vdat)
mm11 = glm(w_male ~ height_c, family="binomial", weights=n*2, data=vdat)
mm12 = glm(w_male ~ flowers_open_c, family="binomial", weights=n*2, data=vdat)
mm13 = glm(w_male ~ flower_size_c, family="binomial", weights=n*2, data=vdat)
mm14 = glm(w_male ~ spur_length_c, family="binomial", weights=n*2, data=vdat)

AICtab = AIC(mm0, mm, mm1, mm2, mm3, mm4, mm5, mm6, mm7, mm8, mm9, mm10, mm11, mm12, mm13, mm14)
AICtab[order(AICtab$AIC),]

mmsum = list(mm, signif(r.squaredGLMM(mm)[1,1], 3))

Gotland_mm = save(mmsum, file="Gotland_mm.RData")

# Fruit set analysis
wdat = na.omit(subset(dat, select = c(n, visited, w_female, fruit_number, flowers)))
wdat = wdat[wdat$fruit_number>0,]

summary(wdat)
hist(wdat$fruit_number)
hist(wdat$fruit_number/wdat$flowers)

mw0 = glm(fruit_number/flowers ~ 1, "binomial", weights=flowers, dat=wdat)
mw = glm(fruit_number/flowers ~ w_female, "binomial", weights=flowers, dat=wdat)

AIC(mw0, mw)

summary(mw)

mwsum = list(mw, signif(r.squaredGLMM(mw)[1,1], 3))

Gotland_mw = save(mwsum, file="Gotland_mw.RData")

# Compute selection gradients
par(mfrow=c(1,1))

# Visitation
rr = summary(mv)$coef[,1]
pd = data.frame(Intercept=rep(1, nrow(dat)), subset(dat, select = names(rr[-1])))
predvis = invlogit(as.matrix(pd) %*% as.matrix(rr))
hist(predvis)

# Poll
rr = summary(mf)$coef[,1]
pd = data.frame(Intercept=rep(1, nrow(dat)), subset(dat, select = names(rr[-1])))
predpoll = invlogit(as.matrix(pd) %*% as.matrix(rr))
predpoll = predpoll*predvis
hist(predpoll)

# Fruitset
rr = summary(mw)$coef[,1]
pd = data.frame(Intercept=rep(1, nrow(dat)), predpoll)
predfruitset = invlogit(as.matrix(pd) %*% as.matrix(rr))
hist(predfruitset)

# Fruits
predfruits = predfruitset*dat$flowers
dat$predfruits = predfruits
hist(predfruits)

# Selection model
relfit = predfruits/mean(predfruits, na.rm=T)
ms = lm(relfit ~ flowers_open_c + height_c + flower_size_c + spur_length_c, dat=dat)
betas = summary(ms)$coef[,1]

betas*100

# Bootstrapping

betaList = list()

for(i in 1:1000){
  
  # Visitation
  rr = MASS::mvrnorm(1, mu=summary(mv)$coef[,1], Sigma=vcov(mv))
  pd = data.frame(Intercept=rep(1, nrow(dat)), subset(dat, select = names(rr[-1])))
  predvis = invlogit(as.matrix(pd) %*% as.matrix(rr))
  
  # Poll
  rr = MASS::mvrnorm(1, mu=summary(mf)$coef[,1], Sigma=vcov(mf))
  pd = data.frame(Intercept=rep(1, nrow(dat)), subset(dat, select = names(rr[-1])))
  predpoll = invlogit(as.matrix(pd) %*% as.matrix(rr))
  
  predpoll = predpoll*predvis
  
  # Fruitset
  rr = MASS::mvrnorm(1, mu=summary(mw)$coef[,1], Sigma=vcov(mw))
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
betas = betas*100

ses = apply(betaDat, 2, sd)*100
apply(betaDat, 2, quantile, c(0.025, 0.975))*100

save(betas, file="Gotland_betas.RData")
save(ses, file="Gotland_ses.RData")

# Evaluating AIC support ####

AICff = AIC(mv) + AIC(mf) + AIC(mw)
AIC0 = AIC(mv0) + AIC(mf0) + AIC(mw0)
AICff = AIC(mv) + AIC(mf)
AIC0 = AIC(mv0) + AIC(mf0)

AICff
AIC0
aicval = AIC0-AICff

save(aicval, file="Gotland_aic.RData")

par(mfrow=c(1,1))
sel = which(dat$fruit_number>0)
plot(dat$fruit_number[sel], dat$predfruits[sel])
lines(-10:100, -10:100)
r2 = cor(dat$fruit_number[sel], dat$predfruits[sel], "pairwise")^2 # Overall r^2 for fitness function

save(r2, file="Gotland_r2.Rdata")

# Naive analysis ####
obs_relfit = dat$fruit_number/mean(dat$fruit_number, na.rm=T)

naivemod = lm(obs_relfit~height_c+flowers_open_c+flower_size_c+spur_length_c, na=na.exclude, dat=dat)

save(naivemod, file="Gotland_mnaive.RData")

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


