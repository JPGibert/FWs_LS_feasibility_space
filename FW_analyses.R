setwd("~/Desktop/JP/Papers_in_review_submitted/JP_FW_constrained_space/Analyses")

library("PerformanceAnalytics")
library('RColorBrewer')
library("quantreg")
library("igraph")
library("dplyr")

## Making log axes ticks
log10Tck <- function(side, type){
   lim <- switch(side, 
     x = par('usr')[1:2],
     y = par('usr')[3:4],
     stop("side argument must be 'x' or 'y'"))
   at <- floor(lim[1]) : ceiling(lim[2])
   return(switch(type, 
     minor = outer(0:9, 10^(min(at):max(at))),
     major = 10^at,
     stop("type argument must be 'major' or 'minor'")
   ))
}

##Load and prepares dataset
web <- read.csv("Food_web_metadata.csv")
struct <- read.csv("struct_data.csv")
web <- cbind(web,struct)
head(web)
web$LoverS <- web$L/web$S
web$absLat <- abs(web$Lat)
web$logL <- log10(web$L)
web$logS <- log10(web$S)
web$logC <- log10(web$C)
web$logOmnv <- log10(web$Omnv+1)
web$logLat <- log10(web$absLat)
#Replacing normal Latitude by rescaled latitude
web$absLat <- (web$absLat-mean(web$absLat))/sd(web$absLat)
web$LogMod <- log10((web$Modularity-mean(web$Modularity))/sd(web$Modularity)+2)
web$Type <- as.factor(web$Type)
web$Type2 <- as.factor(web$Type2)

## Modularity analysis
setwd("~/Desktop/JP/Papers_in_review_submitted/JP_FW_constrained_space/Analyses/Food_webs")
FWs <- dir()
FW_2 <- rep("NA",length(FWs))
for(i in 1:length(FWs)){
	FW_2[i] <- paste(sub(pattern = "\\_Adj.txt", replacement = "", x = FWs[i]))
}

modularity <- rep(0,length(FWs))
for(i in 1:length(FWs)){
	mat <- as.matrix(read.table(FWs[i]))
	g <- graph_from_adjacency_matrix(mat, mode="directed")
	wtc <- cluster_walktrap(g)
	modularity[i] <- modularity(wtc)
}
web$Modularity <- modularity

###################################################################################################
#### The WEDGE nature of food webs

head(web)

## Beta estimation
mean(web$Basal, na.rm=TRUE)
max(web$Basal, na.rm=TRUE)
min(web$Basal, na.rm=TRUE)

pdf('~/Desktop/JP/Papers_in_progress/JP_FW_constrained_space/Figures/Figure_1/Fig_1.pdf', useDingbat=FALSE, width=10, height=5)

par(mfrow=c(1,2), mar=c(0.5,1,2,2), oma=c(3,3,2,3))
web2 <- web

# FIGURE 1b
plot(L~S,log="xy",data=web2,pch=c(21,22,23,24)[web2$Type], col="black",bg=brewer.pal(n = 8, name = "Set1")[2:5][web2$Type], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(6,300), cex=1.2)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
 	abline(a=0,b=1.7, lwd=2, lty=2)
 	abline(a=0,b=1, lwd=2)
 	abline(a=0,b=2, lwd=2)

## FIGURE 1c
qreg <- rq(log(L)~log(S),data=web,tau=seq(0.1,0.90, by=0.05))
plot(L~S,log="xy",data=web2,pch=c(21,22,23,24)[web2$Type], col="black",bg=brewer.pal(n = 8, name = "Set1")[2:5][web2$Type], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(6,300), cex=1.2)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
 	abline(a=0,b=1.7, lwd=2, col="gray")
 	abline(a=0,b=1, lwd=2, col="gray")
 	#abline(a=0,b=2, lwd=2)
	qs <- c(0.1,0.5,0.90)
 for( i in 1:length(qs)){
 	abline(rq(log10(L)~log10(S),data=web,tau=qs[i]), lwd=2)
 }
 points(L~S,log="xy",data=web2,pch=c(21,22,23,24)[web2$Type], col="black",bg=brewer.pal(n = 8, name = "Set1")[2:5][web2$Type], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(6,300), cex=1.2)
dev.off()
summary(qreg)
summary(lm(log10(L)~log10(S),data=web))

## FIGURE 2
# Proportion basal decreases intercept, as predicted by theory…
mod <- lm(logL~logS+logOmnv+LogMod+absLat+Top*Type2+TL+Basal,data=web) # Model without TL and Basal suggests existence of Top*Type2 interaction. Mod was standardized then exponentially trnsformed.
summary(mod) ## Very significant Top*Type interaction

mod <- lm(logL~logS+logOmnv+LogMod+absLat+Top*Type2,data=web)
summary(mod) # We drop insignificant terms without losing much (R^2 goes from 88.6% to 88.3%)

# ATT: Multicollinearity exists (non-significance of factors to be taken with caution) -------------
web_corr <- web %>%
		select(logL,logS,logOmnv,LogMod,absLat,Top,TL,Basal) #Type2, non-numeric so can't be used for correlation
chart.Correlation(web_corr, histogram=FALSE, pch=16)
#---------------------------------------------------------------------------------------------------

pdf('~/Desktop/JP/Papers_in_review_submitted/JP_FW_constrained_space/Figures/Figure_2/Fig_2.pdf', useDingbat=FALSE, width=10, height=6)

size=1.25
par(mfrow=c(2,3), mar=c(0.5,1,2,2), oma=c(3,3,2,3))
## Omnivory
pal = colorRampPalette(c("palegreen", "darkgreen"))
pal(10)[as.numeric(cut(web$logOmnv,breaks = 10))]
TL <- summary(web$logOmnv)
TL <- TL[-c(4,7)]

web2 <- web
plot(L~S,log="xy",data=web2,pch=21, col="black",bg=pal(10)[as.numeric(cut(web2$logOmnv,breaks = 10), na.rm=T)], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(8,270), cex=size)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
 	abline(a=0,b=1.7, lty=2)
 	abline(a=0,b=1, lty=2)
for(i in 1:length(TL)){
	pred_data <- data.frame(logL = seq(min(web2$logL),max(web2$logL),length.out=1000), logS = seq(min(web2$logS),max(web2$logS),length.out=1000), logOmnv = rep(TL[i],length.out=1000),Type2=rep("Aquatic",1000),absLat=rep(mean(web2$absLat),1000), Top=rep(mean(web2$Top, na.rm=T),1000),LogMod=rep(mean(web2$LogMod, na.rm=T)))
	mod_pred <- predict(mod, pred_data)
	pred_plot <- 10^t(t(mod_pred))
	Splot <- 10^t(t(pred_data$logS))
	lines(pred_plot~Splot, col=pal(5)[i], lwd=2)
}

## Modularity
pal = colorRampPalette(brewer.pal(9,"Oranges"))
pal(10)[as.numeric(cut(web2$LogMod,breaks = 10))]
Mod <- summary(web2$LogMod)
Mod <- Mod[-c(4,7)]
plot(L~S,log="xy",data=web2,pch=21, col="black",bg=pal(10)[as.numeric(cut(web2$LogMod,breaks = 10), na.rm=T)], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(8,270), cex=size)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
 	abline(a=0,b=1.7, lty=2)
 	abline(a=0,b=1, lty=2)
for(i in 1:length(Mod)){
	pred_data <- data.frame(logL = seq(min(web2$logL),max(web2$logL),length.out=1000), logS = seq(min(web2$logS),max(web2$logS),length.out=1000), LogMod = rep(Mod[i],length.out=1000),Type2=rep("Aquatic",1000),absLat=rep(mean(web2$absLat),1000), Top=rep(mean(web2$Top, na.rm=T),1000),logOmnv = rep(mean(web2$logOmnv, na.rm=T),length.out=1000))
	mod_pred <- predict(mod, pred_data)
	pred_plot <- 10^t(t(mod_pred))
	Splot <- 10^t(t(pred_data$logS))
	lines(pred_plot~Splot, col=pal(5)[i], lwd=2)
}

## Latitude
pal = colorRampPalette(brewer.pal(9,"Blues"))
Lat <- summary(web2$absLat)
Lat <- Lat[-c(4,7)]
linescol <- rev(pal(8)[8-seq(1,length(Lat))]) # To make color lines more visible
plot(L~S,log="xy",data=web2,pch=21, col="black",bg=pal(10)[as.numeric(cut(web2$absLat,breaks = 10))], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(8,270), cex=size)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
 	#abline(a=0,b=2)
 	abline(a=0,b=1.7, lty=2)
 	abline(a=0,b=1, lty=2)
for(i in 1:length(Lat)){
	pred_data <- data.frame(logL = seq(min(web2$logL),max(web2$logL),length.out=1000), logS = seq(min(web2$logS),max(web2$logS),length.out=1000), logOmnv = rep(mean(web2$logOmnv, na.rm=T),length.out=1000),Type2=rep("Aquatic",1000),absLat=rep(Lat[i],1000),Top=rep(mean(web2$Top, na.rm=T),1000),LogMod=rep(mean(web2$LogMod, na.rm=T)))
	mod_pred <- predict(mod, pred_data)
	pred_plot <- 10^t(t(mod_pred))
	Splot <- 10^t(t(pred_data$logS))
	lines(pred_plot~Splot, col=linescol[i], lwd=2)
}

## Fraction Top (Aquatic)
pal = colorRampPalette(brewer.pal(9,"PuRd"))
Topp <- summary(web2$Top)
Topp <- Topp[-c(4,7)]
Topp[c(2,3,4)] <- Topp[c(2,3,4)]+0.1 # To spread lines a bit more
linescol <- rev(pal(7)[7-seq(1,length(Topp))]) # To make color lines more visible
plot(L[which(web2$Type2=="Aquatic")]~S[which(web2$Type2=="Aquatic")],log="xy",data=web2,pch=21, col="black",bg=pal(10)[as.numeric(cut(web2$Top,breaks = 10), na.rm=T)], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(8,270), cex=size)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=2) # bottom
 	abline(a=0,b=1.7, lty=2)
 	abline(a=0,b=1, lty=2)
for(i in 1:length(Lat)){
	pred_data <- data.frame(logL = seq(min(web2$logL),max(web2$logL),length.out=1000), logS = seq(min(web2$logS),max(web2$logS),length.out=1000), logOmnv = rep(mean(web2$logOmnv, na.rm=T),length.out=1000),Type2=rep("Aquatic",1000),absLat=rep(mean(web2$absLat),1000),Top=rep(Topp[i],1000),LogMod=rep(mean(web2$LogMod, na.rm=T)))
	mod_pred <- predict(mod, pred_data)
	pred_plot <- 10^t(t(mod_pred))
	Splot <- 10^t(t(pred_data$logS))
	lines(pred_plot~Splot, col=linescol[i], lwd=2)
} ## Effect is fully reversed in terrestrial ecosystems

## Fraction Top (Terrestrial)
pal = colorRampPalette(brewer.pal(9,"PuRd"))
Topp <- summary(web2$Top)
Topp <- Topp[-c(4,7)]
Topp[c(2,3,4)] <- Topp[c(2,3,4)]+0.1 # To spread lines a bit more
linescol <- rev(pal(7)[7-seq(1,length(Topp))]) # To make color lines more visible
plot(L[which(web2$Type2=="Terrestrial")]~S[which(web2$Type2=="Terrestrial")],log="xy",data=web2,pch=21, col="black",bg=pal(10)[as.numeric(cut(web2$Top,breaks = 10), na.rm=T)], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(8,270), cex=size)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=2) # bottom
 	abline(a=0,b=1.7, lty=2)
 	abline(a=0,b=1, lty=2)
for(i in 1:length(Lat)){
	pred_data <- data.frame(logL = seq(min(web2$logL),max(web2$logL),length.out=1000), logS = seq(min(web2$logS),max(web2$logS),length.out=1000), logOmnv = rep(mean(web2$logOmnv, na.rm=T),length.out=1000),Type2=rep("Terrestrial",1000),absLat=rep(mean(web2$absLat),1000),Top=rep(Topp[i],1000),LogMod=rep(mean(web2$LogMod, na.rm=T)))
	mod_pred <- predict(mod, pred_data)
	pred_plot <- 10^t(t(mod_pred))
	Splot <- 10^t(t(pred_data$logS))
	lines(pred_plot~Splot, col=linescol[i], lwd=2)
} ## Effect is fully reversed in terrestrial ecosystems

## Ecosystem type
Levels <- levels(web2$Type2)
plot(L~S,log="xy",data=web2,pch=21, col="black",bg=c("purple", "orange")[web2$Type2], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(8,270), cex=size)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
 	#abline(a=0,b=2)
 	abline(a=0,b=1.7, lty=2)
 	abline(a=0,b=1, lty=2)
for(i in 1:length(Levels)){
	pred_data <- data.frame(logL = seq(min(web2$logL),max(web2$logL),length.out=1000), logS = seq(min(web2$logS),max(web2$logS),length.out=1000), logOmnv = rep(mean(web2$logOmnv, na.rm=T),length.out=1000),Type2=rep(Levels[i],1000),absLat=rep(mean(web2$absLat),1000),Top=rep(mean(web2$Top, na.rm=T),1000),LogMod=rep(mean(web2$LogMod, na.rm=T)))
	mod_pred <- predict(mod, pred_data)
	pred_plot <- 10^t(t(mod_pred))
	Splot <- 10^t(t(pred_data$logS))
	lines(pred_plot~Splot, col=c("purple", "orange")[i], lwd=2)
}

dev.off()


#THE END
