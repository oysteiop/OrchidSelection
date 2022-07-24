######################################################################################
##### Analysis of trait-performance-fitness relationships in Orchids (2020-2021) #####
######################################################################################

# Summary of analyses

rm(list=ls())

library(reshape2)
library(plyr)
library(glmmTMB)
library(fBasics)

invlogit = function(x) 1/(1+exp(-x))

# load results
load(file="results/betamat.RData")
load(file="results/Gotland_betas.RData")
betamat = as.data.frame(betamat)
betamat[6,] = betas
rownames(betamat)[6]  = "Dactylorhiza_lapponica_Gotland_2020"

load(file="results/SEmat.RData")
load(file="results/Gotland_ses.RData")
SEmat = as.data.frame(SEmat)
SEmat[6,] = ses
rownames(SEmat)[6]  = "Dactylorhiza_lapponica_Gotland_2020"

load(file="results/sslist.RData")
load(file="results/Gotland_ss.RData")
sslist[[6]] = summary_stats
names(sslist)[6] = "Dactylorhiza_lapponica_Gotland_2020"

load(file="results/prlist.RData")
load(file="results/Gotland_pr.RData")
prlist[[6]] = Gotland_pr
names(prlist)[6] = "Dactylorhiza_lapponica_Gotland_2020"
prDat = rbind.fill(lapply(prlist, function(x) as.data.frame(t(x))))
rownames(prDat) = names(prlist)

load(file="results/pmatlist.RData")
load(file="results/Gotland_pmat.RData")
pmatlist[[6]] = P
names(pmatlist)[6] = "Dactylorhiza_lapponica_Gotland_2020"

load(file="results/mvlist.RData")
load(file="results/Gotland_mv.RData")
mvlist[[6]] = mvsum
names(mvlist)[6] = "Dactylorhiza_lapponica_Gotland_2020"

load(file="results/mflist.RData")
load(file="results/Gotland_mf.RData")
mflist[[6]] = mfsum
names(mflist)[6] = "Dactylorhiza_lapponica_Gotland_2020"

load(file="results/mmlist.RData")
load(file="results/Gotland_mm.RData")
mmlist[[6]] = mmsum
names(mmlist)[6] = "Dactylorhiza_lapponica_Gotland_2020"

load(file="results/mwlist.RData")
load(file="results/Gotland_mw.RData")
mwlist[[6]] = mwsum
names(mwlist)[6] = "Dactylorhiza_lapponica_Gotland_2020"

load(file="results/aiclist.RData")
load(file="results/Gotland_aic.RData")
aiclist[[6]] = aicval
names(aiclist)[6] = "Dactylorhiza_lapponica_Gotland_2020"

load(file="results/r2list.RData")
load(file="results/Gotland_r2.RData")
r2list[[6]] = r2
names(r2list)[6] = "Dactylorhiza_lapponica_Gotland_2020"

load(file="results/mnaivelist.RData")
load(file="results/Gotland_mnaive.RData")
mnaivelist[[6]] = naivemod
#names(mnaivelist)[6] = "Dactylorhiza_lapponica_Gotland_2020"
names(mnaivelist) = names(r2list)

load(file="results/rmflist.RData")
load(file="results/Gotland_rmf.RData")
rmflist[[6]] = rmf
names(rmflist)[6] = "Dactylorhiza_lapponica_Gotland_2020"

# Start of analyses ####

# Deposition-removal correlations
round(unlist(rmflist), 2)
median(unlist(rmflist))
range(unlist(rmflist))

# Compile parameter tables

# Visitation
vmat = rbind.fill(lapply(mvlist[1:5], function(x) as.data.frame(t(summary(x[[1]])$coef$cond[,1]))))
vmat[6,] = rbind.fill(lapply(mvlist[6], function(x) as.data.frame(t(summary(x[[1]])$coef[,1]))))
round(vmat, 3)

vsemat = rbind.fill(lapply(mvlist[1:5], function(x) as.data.frame(t(summary(x[[1]])$coef$cond[,2]))))
vsemat[6,] = rbind.fill(lapply(mvlist[6], function(x) as.data.frame(t(summary(x[[1]])$coef[,2]))))
round(vsemat, 3)

lapply(mvlist, function(x) x[[2]]) #R-squared

# Pollen deposition
fmat = rbind.fill(lapply(mflist[1:5], function(x) as.data.frame(t(summary(x[[1]])$coef$cond[,1]))))
fmat[6,] = rbind.fill(lapply(mflist[6], function(x) as.data.frame(t(summary(x[[1]])$coef[,1]))))
round(fmat, 3)

fsemat = rbind.fill(lapply(mflist[1:5], function(x) as.data.frame(t(summary(x[[1]])$coef$cond[,2]))))
fsemat[6,] = rbind.fill(lapply(mflist[6], function(x) as.data.frame(t(summary(x[[1]])$coef[,2]))))
round(fsemat, 3)

lapply(mflist, function(x) x[[2]]) #R-squared

# Pollinarium removal
mmat = rbind.fill(lapply(mmlist[1:5], function(x) as.data.frame(t(summary(x[[1]])$coef$cond[,1]))))
mmat[6,] = rbind.fill(lapply(mmlist[6], function(x) as.data.frame(t(summary(x[[1]])$coef[,1]))))
round(mmat, 3)

msemat = rbind.fill(lapply(mmlist[1:5], function(x) as.data.frame(t(summary(x[[1]])$coef$cond[,2]))))
msemat[6,] = rbind.fill(lapply(mmlist[6], function(x) as.data.frame(t(summary(x[[1]])$coef[,2]))))
round(msemat, 3)

lapply(mmlist, function(x) x[[2]]) #R-squared

# Fruit set
wmat = rbind.fill(lapply(mwlist[1:5], function(x) as.data.frame(t(summary(x[[1]])$coef$cond[,1]))))
wmat[6,] = rbind.fill(lapply(mwlist[6], function(x) as.data.frame(t(summary(x[[1]])$coef[,1]))))
round(wmat, 3)

wsemat = rbind.fill(lapply(mwlist[1:5], function(x) as.data.frame(t(summary(x[[1]])$coef$cond[,2]))))
wsemat[6,] = rbind.fill(lapply(mwlist[6], function(x) as.data.frame(t(summary(x[[1]])$coef[,2]))))
round(wsemat, 3)

lapply(mwlist, function(x) x[[2]]) #R-squared

# Statistics for fitness function
round(unlist(aiclist),2)
round(unlist(r2list)*100, 1)

# Correlations
pmatlist[[1]]

# Beta summary ####
setable = as.data.frame(round(SEmat[,-1], 5))
svtable = setable^2
meansv = colMeans(svtable)

betatable = as.data.frame(round(betamat[,-1], 1))
colnames(betatable) = c("Flowers", "Height", "Flower size", "Spur length")

means = colMeans(betatable) 
sds = apply(betatable, 2, sd)
sd_corr = sqrt(sds^2 - meansv)

betatable[7,] = means
betatable[8,] = sds
betatable[9,] = sd_corr

rownames(betatable)[7:9] = c("Mean", "SD", "SDcor")

round(betatable, 1)

# Betas plot

cairo_pdf("BetaPlot.pdf", width=10, height=4, fam="Times")
#x11(width=10, height=4)
par(mfrow=c(1,2), mar=c(4,14,2,2))

cols = divPalette(4, "Spectral")

plot(rep(1,6), betatable$Flowers[1:6], pch=16, ylim=c(-100, 100), xlim=c(0.5, 4.5), col=cols[1],
     ylab="Mean-scaled selection gradient (%)", xlab="", las=1,  xaxt="n")
axis(1, at=1:4, labels=F)
text(1:4, par("usr")[3] - 12, srt = 45, adj = 1, cex=.9,
     labels = c("Flowers", "Height", "Flower size", "Spur length"), xpd = TRUE)
points(rep(2,6), betatable$Height[1:6], pch=16, col=cols[2])
points(rep(3,6), betatable$`Flower size`[1:6], pch=16, col=cols[3])
points(rep(4,6), betatable$`Spur length`[1:6], pch=16, col=cols[4])

abline(h=0)


# Pollination reliability plot
#x11(height=4, width=5.5)
par(mar=c(4,4,2,8))
#cols = c("green3", "blue3", "red3", "black")
matplot(prDat[,1], abs(betamat[,2:5]), pch=16, las=1, ylim=c(0,100),
        xlab="",
        ylab="|Mean-scaled selection gradient| (%)",
        col=cols)
mtext("Plants visited (%)", 1, line=2.5)
legend(90, 100, pch=16, col=cols, legend=c("Flowers", "Height", "Flower size", "Spur length"), xpd=T, bty="n")


dev.off()

# Deposition vs. removal ####
fmat
mmat

cor(fmat$flowers_open_c, mmat$flowers_open_c)
cor(fmat$height_c, mmat$height_c)
cor(fmat$flower_size_c, mmat$flower_size_c)
cor(fmat$spur_length_c, mmat$spur_length_c)

cols = divPalette(4, "Spectral")

cairo_pdf("BetaCorFig.pdf", height=4.5, width=4, family="Times")
#x11(height=4.5, width=4)
plot(fmat$flowers_open_c, mmat$flowers_open_c, pch=16, col=cols[1], las=1,
     xlim=c(-5.25,4),
     ylim=c(-2.5,3),
     xlab="",
     ylab="")
mtext("Slope for pollen deposition", 1, line=2.5)
mtext("Slope for pollinarium removal", 2, line=2.5)

segments(fmat$flowers_open_c, mmat$flowers_open_c - msemat$flowers_open_c, 
         fmat$flowers_open_c, mmat$flowers_open_c + msemat$flowers_open_c)
segments(fmat$height_c, mmat$height_c - msemat$height_c, 
         fmat$height_c, mmat$height_c + msemat$height_c)
segments(fmat$flower_size_c, mmat$flower_size_c - msemat$flower_size_c, 
         fmat$flower_size_c, mmat$flower_size_c + msemat$flower_size_c)
segments(fmat$spur_length_c, mmat$spur_length_c - msemat$spur_length_c, 
         fmat$spur_length_c, mmat$spur_length_c + msemat$spur_length_c)

segments(fmat$flowers_open_c - fsemat$flowers_open_c, mmat$flowers_open_c, 
         fmat$flowers_open_c + fsemat$flowers_open_c, mmat$flowers_open_c)
segments(fmat$height_c - fsemat$height_c, mmat$height_c, 
         fmat$height_c + fsemat$height_c, mmat$height_c)
segments(fmat$flower_size_c - fsemat$flower_size_c, mmat$flower_size_c, 
         fmat$flower_size_c + fsemat$flower_size_c, mmat$flower_size_c)
segments(fmat$spur_length_c - fsemat$spur_length_c, mmat$spur_length_c, 
         fmat$spur_length_c + fsemat$spur_length_c, mmat$spur_length_c)

points(fmat$flowers_open_c, mmat$flowers_open_c, pch=16, col=cols[1])
points(fmat$height_c, mmat$height_c, pch=16, col=cols[2])
points(fmat$flower_size_c, mmat$flower_size_c, pch=16, col=cols[3])
points(fmat$spur_length_c, mmat$spur_length_c, pch=16, col=cols[4])

legend("topleft", col=cols, pch=16, cex=.9, bty="n",
       legend=c("Flowers", "Height", "Flower size", "Spur length"))

dev.off()

m = summary(lm(mmat$flowers_open_c~fmat$flowers_open_c))$coef
xx = seq(min(fmat$flowers_open_c), max(fmat$flowers_open_c), 0.01)
yy = m[1,1] + m[2,1]*xx
lines(xx, yy, col=cols[1])

m = summary(lm(mmat$height_c~fmat$height_c))$coef
xx = seq(min(fmat$height_c), max(fmat$height_c), 0.01)
yy = m[1,1] + m[2,1]*xx
lines(xx, yy, col=cols[2])

m = summary(lm(mmat$flower_size_c~fmat$flower_size_c))$coef
xx = seq(min(fmat$flower_size_c), max(fmat$flower_size_c), 0.01)
yy = m[1,1] + m[2,1]*xx
lines(xx, yy, col=cols[3])

m = summary(lm(mmat$spur_length_c~fmat$spur_length_c))$coef
xx = seq(min(fmat$spur_length_c), max(fmat$spur_length_c), 0.01)
yy = m[1,1] + m[2,1]*xx
lines(xx, yy, col=cols[4])

# Misc plots ####
vmat
xx = seq(-1, 1, length.out=10)
cols = 1:6

x11(height=3.5, width=9)
par(mfrow=c(1,3))

plot(xx, seq(0, 1, length.out=10), col="white")
for(i in 1:6){
  yy = invlogit(vmat[i,1]+xx*vmat[i,2])
  lines(xx, yy, col=cols[i])
}

plot(xx, seq(0, 1, length.out=10), col="white")
for(i in 1:6){
  yy = invlogit(vmat[i,1]+xx*vmat[i,3])
  lines(xx, yy, col=cols[i])
}

plot(xx, seq(0, 1, length.out=10), col="white")
for(i in 1:6){
  yy = invlogit(vmat[i,1]+xx*vmat[i,4])
  lines(xx, yy, col=cols[i])
}

x11(height=7, width=7)
par(mfrow=c(2,2))

plot(xx, seq(0, 1, length.out=10), col="white")
for(i in 1:6){
  yy = invlogit(fmat[i,1]+xx*fmat[i,2])
  lines(xx, yy, col=cols[i])
}

plot(xx, seq(0, 1, length.out=10), col="white")
for(i in 1:6){
  yy = invlogit(fmat[i,1]+xx*fmat[i,3])
  lines(xx, yy, col=cols[i])
}

plot(xx, seq(0, 1, length.out=10), col="white")
for(i in 1:6){
  yy = invlogit(fmat[i,1]+xx*fmat[i,4])
  lines(xx, yy, col=cols[i])
}

plot(xx, seq(0, 1, length.out=10), col="white")
for(i in 1:6){
  yy = invlogit(fmat[i,1]+xx*fmat[i,5])
  lines(xx, yy, col=cols[i])
}

# Beta naive summary

naivemat = data.frame(
  summary(mnaivelist[[1]])$coef[-1,1],
  summary(mnaivelist[[2]])$coef[-1,1],
  summary(mnaivelist[[3]])$coef[-1,1],
  summary(mnaivelist[[6]])$coef[-1,1]
)
colnames(naivemat) = names(mnaivelist)[c(1:3, 6)]

naivemat = as.data.frame(t(naivemat))
naivemat

betatable

cor(betatable$Flowers[c(1:3,6)], naivemat$flowers_open_c*100)
cor(betatable$Height[c(1:3,6)], naivemat$height_c*100)
cor(betatable$`Flower size`[c(1:3,6)], naivemat$flower_size_c*100)
cor(betatable$`Spur length`[c(1:3,6)], naivemat$spur_length_c*100)

x11()
plot(betatable$`Spur length`[c(1:3,6)], naivemat$spur_length_c*100,
     xlim=c(-50, 100), ylim=c(-50, 200))
lines(-100:200, -100:200)
points(betatable$`Flower size`[c(1:3,6)], naivemat$flower_size_c*100, col="blue", pch=16)
points(betatable$Height[c(1:3,6)], naivemat$height_c*100, col="brown", pch=16)
points(betatable$Flowers[c(1:3,6)], naivemat$flowers_open_c*100, col="green", pch=16)
abline(h=0, lty=2)
abline(v=0, lty=2)


