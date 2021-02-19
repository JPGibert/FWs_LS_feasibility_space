### clear all data
rm(list=ls())


setwd("~/Desktop/JP/Papers_in_review_submitted/JP_FW_constrained_space/Manuscript/x_Biology_Letters/2_Resubmission")

library("PerformanceAnalytics")
library('RColorBrewer')
library("quantreg")
library("igraph")
library("car")
library("yhat")
library("tidyverse")
library("mgcv")

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
struct <- read.csv("struct_data_New.csv")
web <- cbind(web,struct)
web<-filter(web, !FW.name=="Las Cuevas") %>%
  mutate(Links=ifelse(FW.name=="Carpinter\x92a ", L, Links), Taxa=ifelse(FW.name=="Carpinter\x92a ", S, Taxa))
web$S<-web$Taxa 
web$L<-web$Links 
web$C<-web$L/(web$S^2)


###########################
## Modularity analysis

# The path here should link to the folder where all food webs text files are present.
setwd("~/Desktop/JP/Papers_in_review_submitted/JP_FW_constrained_space/Manuscript/x_Biology_Letters/2_Resubmission/Food_webs")
FWs <- dir(pattern="*.txt") #DW take only .txt files (i.e., omit README)
FWs <- FWs[-36] #DW remove "LasCuevas_Adj.txt"
FW_2 <- rep("NA",length(FWs))
for(i in 1:length(FWs)){
  FW_2[i] <- paste(sub(pattern = "\\_Adj.txt", replacement = "", x = FWs[i]))
}

FWs <- FWs[-1]
modularity <- rep(0,length(FWs))
for(i in 1:length(FWs)){
  mat <- as.matrix(read.table(FWs[i]))
  g <- graph_from_adjacency_matrix(mat, mode="directed")
  wtc <- cluster_walktrap(g)
  modularity[i] <- modularity(wtc)
}
web$Modularity <- modularity
web$LoverS <- web$L/web$S
web$absLat <- abs(web$Lat)
###########################


#DW variable transformations
web$logL <- log10(web$L)
web$logS <- log10(web$S)
web$logC <- log10(web$C)
web$logOmnv <- log10(web$Omnv+1)
web$logTop <- log10(web$Top+1)
web$logTL <- log10(web$TL)
web$logMod <- log10(web$Modularity+1)
web$Type <- as.factor(web$Type)
web$Type2 <- as.factor(web$Type2)



#---------------------------------------------------------------------------------------------------
# FIGURE 1
#---------------------------------------------------------------------------------------------------

sd(web2$Basal, na.rm=TRUE)/sqrt(length(web2$Basal)-1)

web2 <- web
beta_quantiles<-quantile(web2$Basal, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=T) #DW added beta quantiles

#pdf('~/Desktop/JP/Papers_in_review_submitted/JP_FW_constrained_space/Figures/Figure_1/New/Fig_1.pdf', useDingbat=FALSE, width=10, height=5)

par(mfrow=c(1,2), mar=c(0.5,1,2,2), oma=c(3,3,2,3))
# FIGURE 1b
color_palette<-rev(brewer.pal(n = 6, name = "YlOrRd")[-1]) #DW added color palette for beta quantiles
plot(L~S,log="xy",data=web2,pch=c(21,22,23,24)[web2$Type], col="black",bg=brewer.pal(n = 8, name = "Set1")[2:5][web2$Type], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(6,300), cex=1.2)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
 	# abline(a=log10((1-0.73)),b=2, lwd=2, lty=2)
 	for( i in 1:length(beta_quantiles)){
 	  abline(a=log10((1-beta_quantiles[i])), b=2, lwd=2, col=rev(color_palette)[i])
 	}
 	#abline(a=0,b=1, col="grey", lwd=2)
 	abline(a=0,b=2, lwd=2)
 	curve(x-1, 1.001, 1000, add = TRUE, type = "l",ylab = NULL, log = FALSE, n=5000, lwd=2)

## FIGURE 1c
color_palette_2<-brewer.pal(n = 8, name = "YlGnBu")[c(3,5,8)] #DW added color palette for regression quantiles
color_palette_2<-colorRampPalette(color_palette_2)(4)
	qs <- c(0.1,0.9)
qreg <- rq(log10(L)~log10(S),data=web2,tau=qs) #DW changed 'web' to 'web2' | changed 'seq(0.1,0.90, by=0.05)' to 'c(0.1, 0.5, 0.9)' | changed 'log' to 'log10'
plot(L~S,log="xy",data=web2,pch=c(21,22,23,24)[web2$Type], col="black",bg=brewer.pal(n = 8, name = "Set1")[2:5][web2$Type], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(6,300), cex=1.2)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
 	abline(a=log10((1-beta_quantiles[5])), b=2, lwd=2, col="gray") #DW changed to intercept version with 95% beta quantile
 	#abline(a=0,b=1, lwd=2, col="gray")
 for( i in 1:length(qs)){
 	abline(rq(log10(L)~log10(S),data=web2,tau=qs[i]), lwd=3, col=color_palette_2[i]) #DW changed 'web' to 'web2'
 }

curve(x-1, 1.001, 1000, add = TRUE, type = "l",ylab = NULL, log = FALSE, n=5000, lwd=2, col="grey")
 points(L~S,data=web2,pch=c(21,22,23,24)[web2$Type], col="black",bg=brewer.pal(n = 8, name = "Set1")[2:5][web2$Type], xlab="", ylab="",ylim=c(10,4000),xlim=c(6,300), cex=1.2) #DW removed 'axes=FALSE' and 'log="xy"'

dev.off()

summary(qreg,se = "nid")
summary(lm(log10(L)~log10(S),data=web2)) #DW changed 'web' to 'web2'

#---------------------------------------------------------------------------------------------------
# MODEL SELECTION in MuMIn
#---------------------------------------------------------------------------------------------------
web2 <- web[-18,] # dredge cannot deal with NAs so rows containing them are removed

library("MuMIn")
# Dredge shows that Basal and TL can be dropped as both less important and models with lowest AIC don't include them.
full.model <- lm(logL~logS+logOmnv+logMod+absLat+logTop*Type2+TL+Basal,data=web2, na.action = "na.fail")       
dd <- dredge(full.model)
mods <- get.models(dd, subset=TRUE)
get.models(dd, subset = delta < 2)
imp_dd <- importance(dd)
barplot(t(imp_dd), main="Resp")

# Best model contains all variables in model presented in main text (deltaAICc>2 with second best )
full.model <- lm(logL~logS+logOmnv+logMod+absLat+logTop*Type2,data=web2, na.action = "na.fail")       
dd <- dredge(full.model)
imp_dd <- importance(dd)

#pdf('~/Desktop/JP/Papers_in_review_submitted/JP_FW_constrained_space/Manuscript/x_Biology_Letters/2_Resubmission/V_1/Cov_figs/Fig_3.pdf', useDingbat=FALSE, width=10, height=6)
barplot(t(imp_dd), main="", ylab="Importance")
#dev.off()

#---------------------------------------------------------------------------------------------------
# COMMONALITY ANALYSES (How much multicollinearity is actualy there?)
#---------------------------------------------------------------------------------------------------

# ATT: Multicollinearity exists (non-significance of factors to be taken with caution) -------------
web_corr <- web %>%
		select(logL,logS,logOmnv,logMod,absLat,Top,TL,Basal) #Type2, non-numeric so can't be used for correlation
chart.Correlation(web_corr, histogram=FALSE, pch=16)


#DW Commonality analysis -----

#DW standardize predictor variables for commonality analysis
standardized_web<-web %>%
  mutate(across(c("Top", "TL", "Basal", "logS", "absLat", "logOmnv", "logTop", "logTL", "logMod"), function(val){(val-mean(val, na.rm=T))/sd(val, na.rm=T)}))

#DW regressions for both unstandardized and standardized predictors
S_mod<-lm(logL~logS, data=standardized_web) #DW model with only species richness (S) as a predictor
summary(S_mod)

stan_mod<-lm(logL~logS+logOmnv+logMod+absLat+logTop*Type2, data=standardized_web) #DW main model with standardized predictor variables
stan_full_mod<-lm(logL~logS+logOmnv+logMod+absLat+logTop*Type2+logTL+Basal, data=standardized_web) #DW full model with standardized predictor variables

#DW Commonality analysis
reg_model<-'stan_mod'  #DW choose model to analyze

vif_df<-data.frame(VIF=vif(eval(parse(text=reg_model)))) %>% 
  tibble::rownames_to_column(var="Variable")
CA<-regr(eval(parse(text=reg_model))) 
mult_R2<-CA$LM_Output$r.squared
adj_R2<-CA$LM_Output$adj.r.squared
CA_summary<-data.frame(CA$Commonality_Data$CCTotalbyVar) %>%
  tibble::rownames_to_column(var="Variable") %>%
  mutate(U_percent=Unique/mult_R2*100) #DW This provides the % variance explained by the UNIQUE contributions of each variable, not the COMMON contributions which include joint contributions with all other variables
f<-CA$LM_Output$fstatistic
p<-pf(f[1],f[2],f[3],lower.tail=F)

CA$LM_Output #DW Multiple regression summary

summary_table<-data.frame(CA$LM_Output$coef) %>%
  tibble::rownames_to_column(var="Variable") %>%
  mutate(Variable=str_remove(Variable, "Terrestrial"), mult_R2=mult_R2, adj_R2=adj_R2, model_p_val=p) %>%
  left_join(vif_df) %>%
  left_join(CA_summary) %>%
  mutate(across(Estimate:U_percent, ~round(.x, 3))) %>%
  dplyr::select(mult_R2, adj_R2, model_p_val, Variable, Estimate, Std..Error, Pr...t.., VIF, Unique, U_percent, Common, Total)
summary_table #DW Multiple regression summary with VIFs and Commonality Coefficients

# write.csv(summary_table, "/Users/djw68/Desktop/Academic/postdoc/postdoc_duke/duke_work/projects/2020/Food_webs/FWs_LS_feasibility_space-master/stats_FullMod_summary.csv") #DW export regression table


#---------------------------------------------------------------------------------------------------
# FIGURE 2
#---------------------------------------------------------------------------------------------------

#pdf('~/Desktop/JP/Papers_in_review_submitted/JP_FW_constrained_space/Figures/Figure_2/New/Fig_2.pdf', useDingbat=FALSE, width=10, height=6)

mod<-lm(logL~logS+logOmnv+logMod+absLat+logTop*Type2, data=web2) #DW changed to model used above WITHOUT standardized predictors

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
 	curve(x-1, 1.001, 1000, add = TRUE, type = "l",ylab = NULL, log = FALSE, n=5000, lty=2)
 	abline(a=log10((1-beta_quantiles[5])), b=2, lty=2)
 	#abline(a=0,b=1, lty=2)
for(i in 1:length(TL)){
	pred_data <- data.frame(logL = seq(min(web2$logL),max(web2$logL),length.out=1000), logS = seq(min(web2$logS),max(web2$logS),length.out=1000), logOmnv = rep(TL[i],length.out=1000),Type2=rep("Aquatic",1000),absLat=rep(mean(web2$absLat),1000), logTop=rep(mean(web2$logTop, na.rm=T),1000),logMod=rep(mean(web2$logMod, na.rm=T)))
	mod_pred <- predict(mod, pred_data)
	pred_plot <- 10^t(t(mod_pred))
	Splot <- 10^t(t(pred_data$logS))
	lines(pred_plot~Splot, col=pal(5)[i], lwd=2)
}

## Modularity
pal = colorRampPalette(brewer.pal(9,"Oranges"))
pal(10)[as.numeric(cut(web2$logMod,breaks = 10))]
Mod <- summary(web2$logMod)
Mod <- Mod[-c(4,7)]
plot(L~S,log="xy",data=web2,pch=21, col="black",bg=pal(10)[as.numeric(cut(web2$logMod,breaks = 10), na.rm=T)], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(8,270), cex=size)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
	curve(x-1, 1.001, 1000, add = TRUE, type = "l",ylab = NULL, log = FALSE, n=5000, lty=2)
 	abline(a=log10((1-beta_quantiles[5])), b=2, lty=2)
for(i in 1:length(Mod)){
	pred_data <- data.frame(logL = seq(min(web2$logL),max(web2$logL),length.out=1000), logS = seq(min(web2$logS),max(web2$logS),length.out=1000), logMod = rep(Mod[i],length.out=1000),Type2=rep("Aquatic",1000),absLat=rep(mean(web2$absLat),1000), logTop=rep(mean(web2$logTop, na.rm=T),1000),logOmnv = rep(mean(web2$logOmnv, na.rm=T),length.out=1000))
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
	curve(x-1, 1.001, 1000, add = TRUE, type = "l",ylab = NULL, log = FALSE, n=5000, lty=2)
 	abline(a=log10((1-beta_quantiles[5])), b=2, lty=2)
for(i in 1:length(Lat)){
	pred_data <- data.frame(logL = seq(min(web2$logL),max(web2$logL),length.out=1000), logS = seq(min(web2$logS),max(web2$logS),length.out=1000), logOmnv = rep(mean(web2$logOmnv, na.rm=T),length.out=1000),Type2=rep("Aquatic",1000),absLat=rep(Lat[i],1000),logTop=rep(mean(web2$logTop, na.rm=T),1000),logMod=rep(mean(web2$logMod, na.rm=T)))
	mod_pred <- predict(mod, pred_data)
	pred_plot <- 10^t(t(mod_pred))
	Splot <- 10^t(t(pred_data$logS))
	lines(pred_plot~Splot, col=linescol[i], lwd=2)
}

## Fraction logTop (Aquatic)
pal = colorRampPalette(brewer.pal(9,"PuRd"))
Topp <- summary(web2$logTop)
Topp <- Topp[-c(4,7)]
Topp[c(2,3,4)] <- Topp[c(2,3,4)]+0.1 # To spread lines a bit more
linescol <- rev(pal(7)[7-seq(1,length(Topp))]) # To make color lines more visible
plot(L[which(web2$Type2=="Aquatic")]~S[which(web2$Type2=="Aquatic")],log="xy",data=web2,pch=21, col="black",bg=pal(10)[as.numeric(cut(web2$logTop,breaks = 10), na.rm=T)], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(8,270), cex=size)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=2) # bottom
	curve(x-1, 1.001, 1000, add = TRUE, type = "l",ylab = NULL, log = FALSE, n=5000, lty=2)
 	abline(a=log10((1-beta_quantiles[5])), b=2, lty=2)
for(i in 1:length(Lat)){
	pred_data <- data.frame(logL = seq(min(web2$logL),max(web2$logL),length.out=1000), logS = seq(min(web2$logS),max(web2$logS),length.out=1000), logOmnv = rep(mean(web2$logOmnv, na.rm=T),length.out=1000),Type2=rep("Aquatic",1000),absLat=rep(mean(web2$absLat),1000),logTop=rep(Topp[i],1000),logMod=rep(mean(web2$logMod, na.rm=T)))
	mod_pred <- predict(mod, pred_data)
	pred_plot <- 10^t(t(mod_pred))
	Splot <- 10^t(t(pred_data$logS))
	lines(pred_plot~Splot, col=linescol[i], lwd=2)
} ## Effect is fully reversed in terrestrial ecosystems

## Fraction logTop (Terrestrial)
pal = colorRampPalette(brewer.pal(9,"PuRd"))
Topp <- summary(web2$logTop)
Topp <- Topp[-c(4,7)]
Topp[c(2,3,4)] <- Topp[c(2,3,4)]+0.1 # To spread lines a bit more
linescol <- rev(pal(7)[7-seq(1,length(Topp))]) # To make color lines more visible
plot(L[which(web2$Type2=="Terrestrial")]~S[which(web2$Type2=="Terrestrial")],log="xy",data=web2,pch=21, col="black",bg=pal(10)[as.numeric(cut(web2$logTop,breaks = 10), na.rm=T)], xlab="", ylab="", axes=FALSE,ylim=c(10,4000),xlim=c(8,270), cex=size)
box(lwd=2)
axis(1, at=log10Tck('x','major'), tcl= 0.4, lwd=1.5, mgp=c(3,0.10, 0)) # bottom
axis(1, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=1.5) # bottom
axis(2, at=log10Tck('x','major'), tcl= 0.4, las=T, lwd=1.5, mgp=c(3,0.15, 0)) # bottom
axis(2, at=log10Tck('x','minor'), tcl= 0.2, labels=NA, lwd=2) # bottom
	curve(x-1, 1.001, 1000, add = TRUE, type = "l",ylab = NULL, log = FALSE, n=5000, lty=2)
 	abline(a=log10((1-beta_quantiles[5])), b=2, lty=2)
for(i in 1:length(Lat)){
	pred_data <- data.frame(logL = seq(min(web2$logL),max(web2$logL),length.out=1000), logS = seq(min(web2$logS),max(web2$logS),length.out=1000), logOmnv = rep(mean(web2$logOmnv, na.rm=T),length.out=1000),Type2=rep("Terrestrial",1000),absLat=rep(mean(web2$absLat),1000),logTop=rep(Topp[i],1000),logMod=rep(mean(web2$logMod, na.rm=T)))
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
	curve(x-1, 1.001, 1000, add = TRUE, type = "l",ylab = NULL, log = FALSE, n=5000, lty=2)
 	abline(a=log10((1-beta_quantiles[5])), b=2, lty=2)
for(i in 1:length(Levels)){
	pred_data <- data.frame(logL = seq(min(web2$logL),max(web2$logL),length.out=1000), logS = seq(min(web2$logS),max(web2$logS),length.out=1000), logOmnv = rep(mean(web2$logOmnv, na.rm=T),length.out=1000),Type2=rep(Levels[i],1000),absLat=rep(mean(web2$absLat),1000),logTop=rep(mean(web2$logTop, na.rm=T),1000),logMod=rep(mean(web2$logMod, na.rm=T)))
	mod_pred <- predict(mod, pred_data)
	pred_plot <- 10^t(t(mod_pred))
	Splot <- 10^t(t(pred_data$logS))
	lines(pred_plot~Splot, col=c("purple", "orange")[i], lwd=2)
}

dev.off()


#---------------------------------------------------------------------------------------------------
# REMOVAL OF TOP PREDATORS
#---------------------------------------------------------------------------------------------------
## Define function to find props basal, intermediate, top
## Normal food webs
TL_calculator <- function(MAT){
	# Find dimensions of matrix
	dims <- dim(MAT)
	# Crate array of trophic levels and basal, intermediate and top, as well as generality and vulnerability
	gen_vul <- rep(0,0)
	props <- rep(0,0,0)
	TLs <- rep(0,dims[1])
	# Find primary producers, assign TL=1 and add to basal tally
	TLs[rowSums(MAT)==0] <- 1
	props[1] <- sum(TLs)
	# Now calculate all other TLs
	for(j in 1:100){ 
		# We need to iterate several times because the first time the TLs array will only have 0s and 1s
		# and as it starts filling up the measures of TL will change.
		for(i in seq(1:dims[1])[TLs!=1]){
			TLs[i] = 1 + sum(MAT[i,]%*%t(t(TLs)) * (1/sum(MAT[i,]))) 
		} 
	}
	# Now calculate top and intermediate
	# Top predators are not eaten (so marginal column sum == 0), and can't be herbivores (thus TLs>2)
	top_pred <- which(colSums(MAT)[TLs>2]==0)
	props[3] <- sum(colSums(MAT)[TLs>2]==0) 
	props[2] <- dims[1]-props[1]-props[3]
	# Calculate average generality and vulnerability (here a standardized measure of variation in gen and vul)
	gen_vul[1] <- sd(colSums(MAT))/mean(colSums(MAT))
	gen_vul[2] <- sd(rowSums(MAT))/mean(rowSums(MAT))
	# Return TLs
	return(list(TLs, props/dims[1], gen_vul, dim(MAT)[1], top_pred))
}

Post_C <- rep(0,length(FWs))
prop_i <- rep(0,length(FWs))
prop_f <- rep(0,length(FWs))
for(i in 2:length(FWs)){
	mat <- t(as.matrix(read.table(FWs[i]))) # Needs to be transposed so that 1s means i eats j.
	desc <- TL_calculator(mat)
  	prop_i[i] <- desc[[2]][3]
  	if(identical(desc[[5]], integer(0))){ # If there are no top species, keep matrix
  		mat_remove <- mat
  	}else{ # If there are top species, remove those
  		# Remove
  		mat_remove <- mat[-desc[[5]],-desc[[5]]]
  	}
  	# Calculate new Connectance
	Post_C[i] <- sum(mat_remove)/(dim(mat_remove)[1]^2)
	desc <- TL_calculator(mat_remove)
	prop_f[i] <- desc[[2]][3]
}


web$Post_C <- Post_C
web$propdif <- prop_f-prop_i


## Figs in Appendix

#pdf("~/Desktop/JP/Papers_in_review_submitted/JP_FW_constrained_space/Manuscript/x_Biology_Letters/2_Resubmission/V_1/Cov_figs/Fig1.pdf")

hist(web$propdif, main="", xlab="Top_final - Top_initial")
#dev.off()

web$dist <- web$Post_C-web$C
dist <- web$dist[-1]
positive <- web[web$propdif>0,] 
negative <- web[web$propdif<0,] 

#pdf("~/Desktop/JP/Papers_in_review_submitted/JP_FW_constrained_space/Manuscript/x_Biology_Letters/2_Resubmission/V_1/Cov_figs/Fig2.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(positive$dist[order(positive$dist)],pch=16, col=c("purple", "orange")[positive$Type2], ylab="Change in Connectance", main="Food webs where prop Top increases")
abline(a=0,b=0)
plot(negative$dist[order(negative$dist)],pch=16, col=c("purple", "orange")[negative$Type2], ylab="Change in Connectance", main="Food webs where prop Top decreases")
abline(a=0,b=0)
#dev.off()


#THE END
