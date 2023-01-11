#==========================================================
# R script: Tree biomass allometric model.
#
# Copyright (C) March 2012:
# Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
# License: GPL 3 (see license.txt)
#
# The following code is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY.
#
# Reference: Vieilledent G., Vaudry R., Andriamanohisoa S. F. D.,
# Rakotonarivo O. S., Randrianasolo H. Z., Razafindrabe H. N., Bidaud
# Rakotoarivony C., Ebeling J. and Rasamoelina M. 2012. A universal
# approach to estimate biomass and carbon stock in tropical forests
# using generic allometric models. Ecological Applications. 22(2):
# 572-583.
#
#==========================================================

source(file="Legend.R")
library(gdata) # drop.levels

# Site characteristics
Site.Char <- as.data.frame(matrix(nrow=5,ncol=2))
names(Site.Char) <- c("SiteName","ForestType")
Site.Char$SiteName <- c("Bealanana","Fort-Dauphin-Dry","Fort-Dauphin-Moist","Ivohibe","Fandriana")
Site.Char$ForestType <- c("moist-wet","dry","moist-wet","moist-wet","moist-wet")
Levels.Site <- Site.Char$SiteName
nsite <- length(Levels.Site)
Levels.FType <- levels(as.factor(Site.Char$ForestType))
nftype <- length(Levels.FType)

for (i in 1:nftype) {
  
##############################################################
#/////////////////////////////////////////////////////////////
##### Estimating dry density, wet density and moisture content

#####################################
# Mean residual moisture after drying
Data.Moisture.0 <- read.table(file="DataMoisture-G.csv",header=TRUE,sep="\t")
Data.Moisture.0$FType[Data.Moisture.0$Site=="Fort-Dauphin-Dry"] <- "dry"
Data.Moisture.0$FType[Data.Moisture.0$Site!="Fort-Dauphin-Dry"] <- "moist-wet"
Data.Moisture <- drop.levels(subset(Data.Moisture.0,
                           Data.Moisture.0$FType==Levels.FType[i]))
Moisture.mean <- mean(Data.Moisture$Moisture) 

###########################
# Correcting the dry weight
Data.Dens.0 <- read.table(file="DataLab-G.csv",header=TRUE,sep="\t")
Data.Dens.0$FType[Data.Dens.0$Site=="Fort-Dauphin-Dry"] <- "dry"
Data.Dens.0$FType[Data.Dens.0$Site!="Fort-Dauphin-Dry"] <- "moist-wet"
Data.Dens <- drop.levels(subset(Data.Dens.0,Data.Dens.0$FType==Levels.FType[i]))
names(Data.Dens)
Data.Dens$ExactDryWeight <- Data.Dens$DryWeight-(Data.Dens$DryWeight*Moisture.mean)/100

###################
# Exact dry density
Data.Dens$DryDensity <- Data.Dens$ExactDryWeight/Data.Dens$MoistVol2
# By tree
MatTreeLab <- as.data.frame(tapply(Data.Dens$DryDensity,
                                   Data.Dens$Tree,mean,na.rm=TRUE)) 
MatTreeLab$Tree <- row.names(MatTreeLab)
names(MatTreeLab) <- c("DryDensityTree","Tree")
# By species
MatSpeciesLab <- as.data.frame(tapply(Data.Dens$DryDensity,
                                   Data.Dens$Species,mean,na.rm=TRUE)) 
MatSpeciesLab$Species <- row.names(MatSpeciesLab)
names(MatSpeciesLab) <- c("DryDensitySpecies","Species")
# By genus
MatGenusLab <- as.data.frame(tapply(Data.Dens$DryDensity,
                                    Data.Dens$Genus,mean,na.rm=TRUE))
MatGenusLab$Genus <- row.names(MatGenusLab)
names(MatGenusLab) <- c("DryDensityGenus","Genus")
# Plot
pdf(eval(parse(text=paste("file=\"./",Levels.FType[i],"/DryDensity-by-Genus-",Levels.FType[i],".pdf\"",sep=""))))
par(cex=1.4)
boxplot(Data.Dens$DryDensity~substr(Data.Dens$Genus,1,3),
        ylim=c(0,1.5),
        main="Dry density",
        xlab="Genus",
        las=3,
        ylab=expression(paste("Wood density (",g%.%cm^{3},")",sep="")))
abline(h=mean(Data.Dens$DryDensity,na.rm=TRUE),col="red")
dev.off()

#############
# Wet density
Data.Dens$WetDensity <- Data.Dens$MoistWeight/Data.Dens$MoistVol1
# ! strange data with Data.Dens$WetDensity > 1.5 !!
Data.Dens$WetDensity[Data.Dens$WetDensity>1.5] <- NA
# By tree
MatTreeLab$WetDensityTree <- tapply(Data.Dens$WetDensity,
                                    Data.Dens$Tree,mean,
                                    na.rm=TRUE)
# By species
MatSpeciesLab$WetDensitySpecies <- tapply(Data.Dens$WetDensity,
                                    Data.Dens$Species,mean,
                                    na.rm=TRUE)
# By genus
MatGenusLab$WetDensityGenus <- tapply(Data.Dens$WetDensity,
                                      Data.Dens$Genus,mean,
                                      na.rm=TRUE)
# Plot
pdf(eval(parse(text=paste("file=\"./",Levels.FType[i],"/WetDensity-by-Genus-",Levels.FType[i],".pdf\"",sep=""))))
par(cex=1.4)
boxplot(Data.Dens$WetDensity~substr(Data.Dens$Genus,1,3),
        ylim=c(0,1.5),
        main="Fresh density",
        xlab="Genus",
        las=3,
        ylab=expression(paste("Wood density (",g%.%cm^{3},")",sep="")))
abline(h=mean(Data.Dens$WetDensity,na.rm=TRUE),col="red")
dev.off()

#boxplot(Data.Dens$WetDensity~substr(Data.Dens$Genus,1,3))
#boxplot(Data.Dens$DryDensity~substr(Data.Dens$Genus,1,3),
#        add=TRUE)

##################
# Moisture content
Data.Dens$Moisture <- 100*(Data.Dens$WetDensity-Data.Dens$DryDensity)/Data.Dens$WetDensity
# ! strange data with Data.Dens$WetDensity <= Data.Dens$DryDensity and Moisture <=0
Data.Dens$Moisture[Data.Dens$Moisture<=0] <- NA
# By tree
MatTreeLab$MoistureTree <- tapply(Data.Dens$Moisture,Data.Dens$Tree,mean,na.rm=TRUE)
# By species
MatSpeciesLab$MoistureSpecies <- tapply(Data.Dens$Moisture,Data.Dens$Species,mean,na.rm=TRUE)
# By genus
MatGenusLab$MoistureGenus <- tapply(Data.Dens$Moisture,Data.Dens$Genus,mean,na.rm=TRUE)
# Plot
pdf(eval(parse(text=paste("file=\"./",Levels.FType[i],"/Moisture-by-Genus-",Levels.FType[i],".pdf\"",sep=""))))
par(cex=1.4)
boxplot(Data.Dens$Moisture~substr(Data.Dens$Genus,1,3),
        ylim=c(0,100),
        main="Moisture",
        xlab="Genus",
        las=3,
        ylab="Moisture (% of the fresh weight)")
abline(h=mean(Data.Dens$Moisture,na.rm=TRUE),col="red")
dev.off()

####################################################################
#///////////////////////////////////////////////////////////////////
##### Estimating dry AGB (dAGB) based on individual moisture content

Data.AGB.0 <- read.table(file="DataField-G.csv",header=TRUE,sep="\t")
# Forest-type
Data.AGB.0$FType[Data.AGB.0$Site=="Fort-Dauphin-Dry"] <- "dry"
Data.AGB.0$FType[Data.AGB.0$Site!="Fort-Dauphin-Dry"] <- "moist-wet"
# Selecting the forest-type and removing outliers
Data.AGB <- drop.levels(subset(Data.AGB.0,Data.AGB.0$FType==Levels.FType[i]&
                               !is.na(Data.AGB.0$AGB)&
                               !(Data.AGB.0$AGB>2000&Data.AGB.0$DBH<20)))
Data.AGB <- merge(x=Data.AGB,y=MatTreeLab,by="Tree",x.all=TRUE)
names(Data.AGB)
Data.AGB$dAGB.Tree <- Data.AGB$AGB-(Data.AGB$AGB*Data.AGB$MoistureTree)/100
# For compartiments
Data.AGB$dWbranches.Tree <- Data.AGB$Wbranches-(Data.AGB$Wbranches*Data.AGB$MoistureTree)/100
Data.AGB$dWtrunc.Tree <- Data.AGB$Wtrunc-(Data.AGB$Wtrunc*Data.AGB$MoistureTree)/100
Data.AGB$dWleaves.Tree <- Data.AGB$Wleaves-(Data.AGB$Wleaves*Data.AGB$MoistureTree)/100

#######
#////// Note: number of trees and DBH range by site
sink(file="ntree-by-site.txt")
print(table(Data.AGB.0$Site))
sink()
sink(file="Drange-by-site.txt")
print(tapply(Data.AGB.0$DBH,Data.AGB.0$Site,range,na.rm=TRUE))
sink()

####################################################################
#///////////////////////////////////////////////////////////////////
##### Estimating dry AGB (dAGB) based on genus moisture content

Data.AGB <- merge(x=Data.AGB,y=MatGenusLab,by="Genus",x.all=TRUE)
names(Data.AGB)
Data.AGB$dAGB.Genus <- Data.AGB$AGB-(Data.AGB$AGB*Data.AGB$MoistureGenus)/100
# For compartiments
Data.AGB$dWbranches.Genus <- Data.AGB$Wbranches-(Data.AGB$Wbranches*Data.AGB$MoistureGenus)/100
Data.AGB$dWtrunc.Genus <- Data.AGB$Wtrunc-(Data.AGB$Wtrunc*Data.AGB$MoistureGenus)/100
Data.AGB$dWleaves.Genus <- Data.AGB$Wleaves-(Data.AGB$Wleaves*Data.AGB$MoistureGenus)/100

# Check for difference between dAGB.Tree and dAGB.Genus
plot(Data.AGB$dAGB.Tree,Data.AGB$dAGB.Genus) # this is the same

# We choose dAGB.Tree if it exists and dAGB.Genus else
Data.AGB$dAGB <- Data.AGB$dAGB.Tree
Data.AGB$dAGB[is.na(Data.AGB$dAGB.Tree)] <- Data.AGB$dAGB.Genus[is.na(Data.AGB$dAGB.Tree)]
# For compartiments
Data.AGB$dWbranches <- Data.AGB$dWbranches.Tree
Data.AGB$dWbranches[is.na(Data.AGB$dWbranches.Tree)] <- Data.AGB$dWbranches.Genus[is.na(Data.AGB$dWbranches.Tree)]
Data.AGB$dWtrunc <- Data.AGB$dWtrunc.Tree
Data.AGB$dWtrunc[is.na(Data.AGB$dWtrunc.Tree)] <- Data.AGB$dWtrunc.Genus[is.na(Data.AGB$dWtrunc.Tree)]
Data.AGB$dWleaves <- Data.AGB$dWleaves.Tree
Data.AGB$dWleaves[is.na(Data.AGB$dWleaves.Tree)] <- Data.AGB$dWleaves.Genus[is.na(Data.AGB$dWleaves.Tree)]

###########################################################
#//////////////////////////////////////////////////////////
##### Allometric data

# Log transformation
ldAGB <- Data.AGB$ldAGB <- log(Data.AGB$dAGB) 
lD <- Data.AGB$lD <- log(Data.AGB$DBH)
lH <- Data.AGB$lH <- log(Data.AGB$Htot)
lrho <- Data.AGB$lrho <- log(Data.AGB$DryDensityTree) # Must be DryDensityTree if available !!!
lrho[is.na(lrho)] <- Data.AGB$lrho[is.na(Data.AGB$lrho)] <- log(Data.AGB$DryDensityGenus[is.na(Data.AGB$DryDensityTree)])

# Export raw data
RawData <- data.frame(1:length(ldAGB),Data.AGB$Site,Data.AGB$Species,Data.AGB$Family,exp(lD),exp(lH),
                      Data.AGB$dWbranches,Data.AGB$dWtrunc,Data.AGB$dWleaves,
                      exp(ldAGB),exp(lrho),rep(Levels.FType[i],length(ldAGB)))
names(RawData) <- c("Tree","Locality","Latin_Binomial","Family","Trunk_diameter_cm","Total_height_m",
                    "Dry_weight_branch_kg","Dry_weight_trunk_kg","Dry_weight_leaf_kg",
                    "Dry_total_AGB_kg","Wood_specific_gravity","Forest_type")
RawData <- RawData[order(RawData$Latin_Binomial,RawData$Trunk_diameter_cm),]
RawData$Tree <- c(1:length(ldAGB)) 
write.table(RawData,
            file=eval(parse(text=paste("file=\"./",Levels.FType[i],"/Raw-Data-",Levels.FType[i],".txt\"",sep=""))),
            col.names=TRUE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

########################################################################################
#======================================================================================#
#================================= Model comparison ===================================#

###################################################
#========= Data-frame to record results ==========#
Result.Mat <- as.data.frame(matrix(ncol=11,nrow=9))
names(Result.Mat) <- c("Model","beta0","beta1","beta2",
                       "beta3","beta4","df","RSE","R2","AIC","Bias")
Result.Mat$Model <- c("Brown","Chave.H","Chave.D","Mada.I.1","Mada.I.1-wHD","Mada.I.2","Mada.I.2-wHD",
                      "Mada.II.1","Mada.II.2-WEB")

##########################
#========= Brown ========#
if (Levels.FType[i]=="moist-wet") {
  # Moist forest:
  # log(AGB)=-2.289+2.649*log(D)-0.021*(log(D)^2)

  # Predictions from Brown model
  ldAGB.Brown <- -2.289+2.649*lD-0.021*(lD)^2
}

if (Levels.FType[i]=="dry") {
  # Dry forest:
  # log(AGB)=-0.535*log(10)+log(pi/4)+2*log(D)

  # Predictions from Brown model
  ldAGB.Brown <- -0.535*log(10)+log(pi/4)+2*lD
}

# Graph
dev.set(dev.prev())
plot(ldAGB,ldAGB.Brown,xlim=c(0,10),ylim=c(0,10))
lines(seq(0,10),seq(0,10),col="red")

# df
nobs <- length(Data.AGB$dAGB)
if (Levels.FType[i]=="moist-wet") {
  npar <- 3
}
if (Levels.FType[i]=="dry") {
  npar <- 2
}
df <- nobs-npar

# RSE
# Residuals
Res <- ldAGB-ldAGB.Brown
# Distribution of the residuals
hist(Res) # Strong overestimation of biomass with Brown equation
# Bias
Bias.Tree.Brown <- 100*(1/exp(Res)-1)
Bias <- mean(100*(1/exp(Res)-1)) # Biased estimation of Brown model (biais)

# Residuals Standard Error
RSE <- sd(Res)

# R square
SCR <- sum((ldAGB-ldAGB.Brown)^2)
SCT <- sum((ldAGB-mean(ldAGB))^2)
R.square <- 1-SCR/SCT

# AIC
logL <- sum(dnorm(ldAGB,mean=ldAGB.Brown,sd=sd(Res),log=TRUE))
AIC <- -2*logL+2*npar

# Completing the result matrix
if (Levels.FType[i]=="moist-wet") {
  Result.Mat$beta0[Result.Mat$Model=="Brown"] <- -2.289
  Result.Mat$beta1[Result.Mat$Model=="Brown"] <- 2.649
  Result.Mat$beta2[Result.Mat$Model=="Brown"] <- -0.021
  Result.Mat$beta3[Result.Mat$Model=="Brown"] <- "-"
  Result.Mat$beta4[Result.Mat$Model=="Brown"] <- "-"
  Result.Mat$df[Result.Mat$Model=="Brown"] <- df
  Result.Mat$RSE[Result.Mat$Model=="Brown"] <- RSE
  Result.Mat$R2[Result.Mat$Model=="Brown"] <- R.square
  Result.Mat$AIC[Result.Mat$Model=="Brown"] <- AIC
  Result.Mat$Bias[Result.Mat$Model=="Brown"] <- Bias
}
if (Levels.FType[i]=="dry") {
  Result.Mat$beta0[Result.Mat$Model=="Brown"] <- -0.535*log(10)+log(pi/4)
  Result.Mat$beta1[Result.Mat$Model=="Brown"] <- 2
  Result.Mat$beta2[Result.Mat$Model=="Brown"] <- "-"
  Result.Mat$beta3[Result.Mat$Model=="Brown"] <- "-"
  Result.Mat$beta4[Result.Mat$Model=="Brown"] <- "-"
  Result.Mat$df[Result.Mat$Model=="Brown"] <- df
  Result.Mat$RSE[Result.Mat$Model=="Brown"] <- RSE
  Result.Mat$R2[Result.Mat$Model=="Brown"] <- R.square
  Result.Mat$AIC[Result.Mat$Model=="Brown"] <- AIC
  Result.Mat$Bias[Result.Mat$Model=="Brown"] <- Bias
}

############################
#========= Chave.H ========#

if (Levels.FType[i]=="moist-wet") {
  # Moist forest:
  # log(AGB)=-2.977+log(D^2*H*rho)

  # Predictions from Chave.H model
  ldAGB.Chave.H <- -2.977+log(Data.AGB$DBH^2*Data.AGB$Htot*
                                  Data.AGB$DryDensityGenus)
}

if (Levels.FType[i]=="dry") {
  # Moist forest:
  # log(AGB)=-2.187+0.916*log(D^2*H*rho)

  # Predictions from Chave.H model
  ldAGB.Chave.H <- -2.187+0.916*log(Data.AGB$DBH^2*Data.AGB$Htot*
                                  Data.AGB$DryDensityGenus)
}

# Graph
plot(ldAGB,ldAGB.Chave.H,xlim=c(0,10),ylim=c(0,10))
lines(seq(0,10),seq(0,10),col="red")

# df
nobs <- length(Data.AGB$dAGB)
if (Levels.FType[i]=="moist-wet") {
  npar <- 1
}
if (Levels.FType[i]=="dry") {
  npar <- 2
}
df <- nobs-npar

# RSE
# Residuals
Res <- ldAGB-ldAGB.Chave.H
# Distribution of the residuals
hist(Res) # Strong overestimation of biomass with Chave equation
# Residuals Standard Error
RSE <- sd(Res)
# Bias
Bias.Tree.Chave.H <- 100*(1/exp(Res)-1)
Bias <- mean(100*(1/exp(Res)-1)) # Biased estimation of Chave model (biais)

# R square
SCR <- sum((ldAGB-ldAGB.Chave.H)^2)
SCT <- sum((ldAGB-mean(ldAGB))^2)
R.square <- 1-SCR/SCT

# AIC
logL <- sum(dnorm(ldAGB,mean=ldAGB.Chave.H,sd=sd(Res),log=TRUE))
AIC <- -2*logL+2*npar

# Completing the result matrix
if (Levels.FType[i]=="moist-wet") {
  Result.Mat$beta0[Result.Mat$Model=="Chave.H"] <- -2.977
  Result.Mat$beta1[Result.Mat$Model=="Chave.H"] <- "-"
  Result.Mat$beta2[Result.Mat$Model=="Chave.H"] <- "-"
  Result.Mat$beta3[Result.Mat$Model=="Chave.H"] <- "-"
  Result.Mat$beta4[Result.Mat$Model=="Chave.H"] <- "-"
  Result.Mat$df[Result.Mat$Model=="Chave.H"] <- df
  Result.Mat$RSE[Result.Mat$Model=="Chave.H"] <- RSE
  Result.Mat$R2[Result.Mat$Model=="Chave.H"] <- R.square
  Result.Mat$AIC[Result.Mat$Model=="Chave.H"] <- AIC
  Result.Mat$Bias[Result.Mat$Model=="Chave.H"] <- Bias
}

if (Levels.FType[i]=="dry") {
  Result.Mat$beta0[Result.Mat$Model=="Chave.H"] <- -2.187
  Result.Mat$beta1[Result.Mat$Model=="Chave.H"] <- 0.916
  Result.Mat$beta2[Result.Mat$Model=="Chave.H"] <- "-"
  Result.Mat$beta3[Result.Mat$Model=="Chave.H"] <- "-"
  Result.Mat$beta4[Result.Mat$Model=="Chave.H"] <- "-"
  Result.Mat$df[Result.Mat$Model=="Chave.H"] <- df
  Result.Mat$RSE[Result.Mat$Model=="Chave.H"] <- RSE
  Result.Mat$R2[Result.Mat$Model=="Chave.H"] <- R.square
  Result.Mat$AIC[Result.Mat$Model=="Chave.H"] <- AIC
  Result.Mat$Bias[Result.Mat$Model=="Chave.H"] <- Bias
}

##############################
#========= Chave.D ========#

if (Levels.FType[i]=="moist-wet") {
  # Moist forest:
  # log(AGB)=-1.499+2.148*log(D)+0.207*(log(D))^2-0.0281*(log(D))^3+log(rho)

  # Predictions from Chave.D model
  ldAGB.Chave.D <- -1.499+2.148*lD+0.207*(lD)^2-0.0281*(lD)^3+lrho
}

if (Levels.FType[i]=="dry") {
  # Dry forest:
  # log(AGB)=-0.667+1.784*log(D)+0.207*(log(D))^2-0.0281*(log(D))^3+log(rho)

  # Predictions from Chave.D model
  ldAGB.Chave.D <- -0.667+1.784*lD+0.207*(lD)^2-0.0281*(lD)^3+lrho
}

# Graph
plot(ldAGB,ldAGB.Chave.D,xlim=c(0,10),ylim=c(0,10))
lines(seq(0,10),seq(0,10),col="red")

# df
nobs <- length(Data.AGB$dAGB)
if (Levels.FType[i]=="moist-wet") {
  npar <- 4
}
if (Levels.FType[i]=="dry") {
  npar <- 4
}
df <- nobs-npar

# RSE
# Residuals
Res <- ldAGB-ldAGB.Chave.D
# Distribution of the residuals
hist(Res) # Strong overestimation of biomass with Chave equation
# Residuals Standard Error
RSE <- sd(Res)
# Bias
Bias.Tree.Chave.D <- 100*(1/exp(Res)-1)
Bias <- mean(100*(1/exp(Res)-1)) # Biased estimation of Chave model (biais)

# R square
SCR <- sum((ldAGB-ldAGB.Chave.D)^2)
SCT <- sum((ldAGB-mean(ldAGB))^2)
R.square <- 1-SCR/SCT

# AIC
logL <- sum(dnorm(ldAGB,mean=ldAGB.Chave.D,sd=sd(Res),log=TRUE))
AIC <- -2*logL+2*npar

# Completing the result matrix
if (Levels.FType[i]=="moist-wet") {
  Result.Mat$beta0[Result.Mat$Model=="Chave.D"] <- -1.499
  Result.Mat$beta1[Result.Mat$Model=="Chave.D"] <- 2.148
  Result.Mat$beta2[Result.Mat$Model=="Chave.D"] <- 0.207
  Result.Mat$beta3[Result.Mat$Model=="Chave.D"] <- -0.0281
  Result.Mat$beta4[Result.Mat$Model=="Chave.D"] <- "-"
  Result.Mat$df[Result.Mat$Model=="Chave.D"] <- df
  Result.Mat$RSE[Result.Mat$Model=="Chave.D"] <- RSE
  Result.Mat$R2[Result.Mat$Model=="Chave.D"] <- R.square
  Result.Mat$AIC[Result.Mat$Model=="Chave.D"] <- AIC
  Result.Mat$Bias[Result.Mat$Model=="Chave.D"] <- Bias
}

if (Levels.FType[i]=="dry") {
  Result.Mat$beta0[Result.Mat$Model=="Chave.D"] <- -0.667
  Result.Mat$beta1[Result.Mat$Model=="Chave.D"] <- 1.784
  Result.Mat$beta2[Result.Mat$Model=="Chave.D"] <- 0.207
  Result.Mat$beta3[Result.Mat$Model=="Chave.D"] <- -0.0281
  Result.Mat$beta4[Result.Mat$Model=="Chave.D"] <- "-"
  Result.Mat$df[Result.Mat$Model=="Chave.D"] <- df
  Result.Mat$RSE[Result.Mat$Model=="Chave.D"] <- RSE
  Result.Mat$R2[Result.Mat$Model=="Chave.D"] <- R.square
  Result.Mat$AIC[Result.Mat$Model=="Chave.D"] <- AIC
  Result.Mat$Bias[Result.Mat$Model=="Chave.D"] <- Bias
}


#################################################################
#========= Biomass-diameter-height regression Model I.1 ========#

# Inference
M <- lm(ldAGB~lD+lH+lrho,data=Data.AGB)

# Prediction with Madagascar allometry
ldAGB.Mada.I.1 <- M$fitted

# Graph
plot(ldAGB,as.numeric(ldAGB.Mada.I.1),xlim=c(0,10),ylim=c(0,10))
lines(seq(0,10),seq(0,10),col="red")

# df
nobs <- length(Data.AGB$dAGB)
npar <- 4
df <- nobs-npar

# RSE
# Residuals
Res <- ldAGB-ldAGB.Mada.I.1
# Distribution of the residuals
hist(Res)
# Residuals Standard Error
RSE <- sd(Res)
# Bias
Bias.Tree.Mada.I.1 <- 100*(1/exp(Res)-1)
Bias <- mean(100*(1/exp(Res)-1)) # Biased estimation of our model (biais)

# R square
SCR <- sum((ldAGB-ldAGB.Mada.I.1)^2)
SCT <- sum((ldAGB-mean(ldAGB))^2)
R.square <- 1-SCR/SCT

# AIC
logL <- sum(dnorm(ldAGB,mean=ldAGB.Mada.I.1,sd=sd(Res),log=TRUE))
AIC <- -2*logL+2*npar

# Completing the result matrix
Result.Mat$beta0[Result.Mat$Model=="Mada.I.1"] <- M$coeff[1]+(sd(M$residuals)^2)/2
Result.Mat$beta1[Result.Mat$Model=="Mada.I.1"] <- M$coeff[2]
Result.Mat$beta2[Result.Mat$Model=="Mada.I.1"] <- M$coeff[3]
Result.Mat$beta3[Result.Mat$Model=="Mada.I.1"] <- "-"
Result.Mat$beta4[Result.Mat$Model=="Mada.I.1"] <- M$coeff[4]
Result.Mat$df[Result.Mat$Model=="Mada.I.1"] <- df
Result.Mat$RSE[Result.Mat$Model=="Mada.I.1"] <- RSE
Result.Mat$R2[Result.Mat$Model=="Mada.I.1"] <- R.square
Result.Mat$AIC[Result.Mat$Model=="Mada.I.1"] <- AIC
Result.Mat$Bias[Result.Mat$Model=="Mada.I.1"] <- Bias

# Test
sink(file=eval(parse(text=paste("file=\"./",Levels.FType[i],"/summaryMI1",Levels.FType[i],".txt\"",sep=""))))
print(summary(M))
sink()

#=================#
# With HDallometry

# Construction new data-set (ndt)
if (Levels.FType[i]=="moist-wet") {
  lH <- 1.010+0.547*lD
}
if (Levels.FType[i]=="dry") {
  lH <- log(12.120-(12.120-1.3)*exp(-0.052*Data.AGB$DBH))
}

# Prediction with Madagascar allometry
ldAGB.Mada.I.1.wHD <- predict.lm(M,newdata=as.data.frame(cbind(lD,lH,lrho)))

# Graph
plot(ldAGB,as.numeric(ldAGB.Mada.I.1.wHD),xlim=c(0,10),ylim=c(0,10))
lines(seq(0,10),seq(0,10),col="red")

# df
nobs <- length(Data.AGB$dAGB)
npar <- 4
df <- nobs-npar

# RSE
# Residuals
Res <- ldAGB-ldAGB.Mada.I.1.wHD
# Distribution of the residuals
hist(Res)
# Residuals Standard Error
RSE <- sd(Res)
# Bias
Bias.Tree.Mada.I.1.wHD <- 100*(1/exp(Res)-1)
Bias <- mean(100*(1/exp(Res)-1)) # Biased estimation of our model (biais)

# R square
SCR <- sum((ldAGB-ldAGB.Mada.I.1.wHD)^2)
SCT <- sum((ldAGB-mean(ldAGB))^2)
R.square <- 1-SCR/SCT

# AIC
logL <- sum(dnorm(ldAGB,mean=ldAGB.Mada.I.1.wHD,sd=sd(Res),log=TRUE))
AIC <- -2*logL+2*npar

# Completing the result matrix
if (Levels.FType[i]=="moist-wet") {
  Result.Mat$beta0[Result.Mat$Model=="Mada.I.1-wHD"] <- 1.010
  Result.Mat$beta1[Result.Mat$Model=="Mada.I.1-wHD"] <- 0.547
  Result.Mat$beta2[Result.Mat$Model=="Mada.I.1-wHD"] <- 0.071
  Result.Mat$beta3[Result.Mat$Model=="Mada.I.1-wHD"] <- "-"
  Result.Mat$beta4[Result.Mat$Model=="Mada.I.1-wHD"] <- "-"
  Result.Mat$df[Result.Mat$Model=="Mada.I.1-wHD"] <- df
  Result.Mat$RSE[Result.Mat$Model=="Mada.I.1-wHD"] <- RSE
  Result.Mat$R2[Result.Mat$Model=="Mada.I.1-wHD"] <- R.square
  Result.Mat$AIC[Result.Mat$Model=="Mada.I.1-wHD"] <- AIC
  Result.Mat$Bias[Result.Mat$Model=="Mada.I.1-wHD"] <- Bias
}
if (Levels.FType[i]=="dry") {
  Result.Mat$beta0[Result.Mat$Model=="Mada.I.1-wHD"] <- 12.120
  Result.Mat$beta1[Result.Mat$Model=="Mada.I.1-wHD"] <- 1.300
  Result.Mat$beta2[Result.Mat$Model=="Mada.I.1-wHD"] <- 0.052
  Result.Mat$beta3[Result.Mat$Model=="Mada.I.1-wHD"] <- "-"
  Result.Mat$beta4[Result.Mat$Model=="Mada.I.1-wHD"] <- "-"
  Result.Mat$df[Result.Mat$Model=="Mada.I.1-wHD"] <- df
  Result.Mat$RSE[Result.Mat$Model=="Mada.I.1-wHD"] <- RSE
  Result.Mat$R2[Result.Mat$Model=="Mada.I.1-wHD"] <- R.square
  Result.Mat$AIC[Result.Mat$Model=="Mada.I.1-wHD"] <- AIC
  Result.Mat$Bias[Result.Mat$Model=="Mada.I.1-wHD"] <- Bias
}

#################################################################
#========= Biomass-diameter-height regression Model I.2 ========#

# Inference
M <- lm(ldAGB~I(log(DBH^2*Htot*exp(lrho))),data=Data.AGB)

# Prediction with Madagascar allometry
ldAGB.Mada.I.2 <- M$fitted

# Graph
plot(ldAGB,ldAGB.Mada.I.2,xlim=c(0,10),ylim=c(0,10))
lines(seq(0,10),seq(0,10),col="red")

# df
nobs <- length(Data.AGB$dAGB)
npar <- 2
df <- nobs-npar

# RSE
# Residuals
Res <- ldAGB-ldAGB.Mada.I.2
# Distribution of the residuals
hist(Res)
# Residuals Standard Error
RSE <- sd(Res)
# Bias
Bias.Tree.Mada.I.2 <- 100*(1/exp(Res)-1)
Bias <- mean(100*(1/exp(Res)-1)) # Biased estimation of our model (biais)

# R square
SCR <- sum((ldAGB-ldAGB.Mada.I.2)^2)
SCT <- sum((ldAGB-mean(ldAGB))^2)
R.square <- 1-SCR/SCT

# AIC
logL <- sum(dnorm(ldAGB,mean=ldAGB.Mada.I.2,sd=sd(Res),log=TRUE))
AIC <- -2*logL+2*npar

# Completing the result matrix
Result.Mat$beta0[Result.Mat$Model=="Mada.I.2"] <- M$coeff[1]+(sd(M$residuals)^2)/2
Result.Mat$beta1[Result.Mat$Model=="Mada.I.2"] <- M$coeff[2]
Result.Mat$beta2[Result.Mat$Model=="Mada.I.2"] <- "-"
Result.Mat$beta3[Result.Mat$Model=="Mada.I.2"] <- "-"
Result.Mat$beta4[Result.Mat$Model=="Mada.I.2"] <- "-"
Result.Mat$df[Result.Mat$Model=="Mada.I.2"] <- df
Result.Mat$RSE[Result.Mat$Model=="Mada.I.2"] <- RSE
Result.Mat$R2[Result.Mat$Model=="Mada.I.2"] <- R.square
Result.Mat$AIC[Result.Mat$Model=="Mada.I.2"] <- AIC
Result.Mat$Bias[Result.Mat$Model=="Mada.I.2"] <- Bias

# Test
sink(file=eval(parse(text=paste("file=\"./",Levels.FType[i],"/summaryMI2",Levels.FType[i],".txt\"",sep=""))))
print(summary(M))
sink()

#=================#
# With HDallometry

# Construction new data-set (ndt)
if (Levels.FType[i]=="moist-wet") {
  Htot <- exp(1.132+0.529*lD)
}
if (Levels.FType[i]=="dry") {
  Htot <- 12.12-(12.12-1.3)*exp(-0.052*Data.AGB$DBH)
}

DBH <- Data.AGB$DBH

# Prediction with Madagascar allometry
ldAGB.Mada.I.2.wHD <- predict.lm(M,newdata=as.data.frame(cbind(DBH,Htot,lrho)))

# Graph
plot(ldAGB,as.numeric(ldAGB.Mada.I.2.wHD),xlim=c(0,10),ylim=c(0,10))
lines(seq(0,10),seq(0,10),col="red")

# df
nobs <- length(Data.AGB$dAGB)
npar <- 2
df <- nobs-npar

# RSE
# Residuals
Res <- ldAGB-ldAGB.Mada.I.2.wHD
# Distribution of the residuals
hist(Res)
# Residuals Standard Error
RSE <- sd(Res)
# Bias
Bias.Tree.Mada.I.2.wHD <- 100*(1/exp(Res)-1)
Bias <- mean(100*(1/exp(Res)-1)) # Biased estimation of our model (biais)

# R square
SCR <- sum((ldAGB-ldAGB.Mada.I.2.wHD)^2)
SCT <- sum((ldAGB-mean(ldAGB))^2)
R.square <- 1-SCR/SCT

# AIC
logL <- sum(dnorm(ldAGB,mean=ldAGB.Mada.I.2.wHD,sd=sd(Res),log=TRUE))
AIC <- -2*logL+2*npar

# Completing the result matrix
if (Levels.FType[i]=="moist-wet") {
  Result.Mat$beta0[Result.Mat$Model=="Mada.I.2-wHD"] <- 1.132
  Result.Mat$beta1[Result.Mat$Model=="Mada.I.2-wHD"] <- 0.529
  Result.Mat$beta2[Result.Mat$Model=="Mada.I.2-wHD"] <- "-"
  Result.Mat$beta3[Result.Mat$Model=="Mada.I.2-wHD"] <- "-"
  Result.Mat$beta4[Result.Mat$Model=="Mada.I.2-wHD"] <- "-"
  Result.Mat$df[Result.Mat$Model=="Mada.I.2-wHD"] <- df
  Result.Mat$RSE[Result.Mat$Model=="Mada.I.2-wHD"] <- RSE
  Result.Mat$R2[Result.Mat$Model=="Mada.I.2-wHD"] <- R.square
  Result.Mat$AIC[Result.Mat$Model=="Mada.I.2-wHD"] <- AIC
  Result.Mat$Bias[Result.Mat$Model=="Mada.I.2-wHD"] <- Bias
}
if (Levels.FType[i]=="dry") {
  Result.Mat$beta0[Result.Mat$Model=="Mada.I.2-wHD"] <- 12.12
  Result.Mat$beta1[Result.Mat$Model=="Mada.I.2-wHD"] <- 1.300
  Result.Mat$beta2[Result.Mat$Model=="Mada.I.2-wHD"] <- 0.052
  Result.Mat$beta3[Result.Mat$Model=="Mada.I.2-wHD"] <- "-"
  Result.Mat$beta4[Result.Mat$Model=="Mada.I.2-wHD"] <- "-"
  Result.Mat$df[Result.Mat$Model=="Mada.I.2-wHD"] <- df
  Result.Mat$RSE[Result.Mat$Model=="Mada.I.2-wHD"] <- RSE
  Result.Mat$R2[Result.Mat$Model=="Mada.I.2-wHD"] <- R.square
  Result.Mat$AIC[Result.Mat$Model=="Mada.I.2-wHD"] <- AIC
  Result.Mat$Bias[Result.Mat$Model=="Mada.I.2-wHD"] <- Bias
}

###########################################################
#========= Biomass-diameter regression Model II.1 ========#

# Inference
M <- lm(ldAGB~lD+lrho,data=Data.AGB)

# Prediction with Madagascar allometry
ldAGB.Mada.II.1 <- M$fitted

# Graph
plot(ldAGB,ldAGB.Mada.II.1,xlim=c(0,10),ylim=c(0,10))
lines(seq(0,10),seq(0,10),col="red")

# df
nobs <- length(Data.AGB$dAGB)
npar <- 3
df <- nobs-npar

# RSE
# Residuals
Res <- ldAGB-ldAGB.Mada.II.1
# Distribution of the residuals
hist(Res)
# Residuals Standard Error
RSE <- sd(Res)
# Bias
Bias.Tree.Mada.II.1 <- 100*(1/exp(Res)-1)
Bias <- mean(100*(1/exp(Res)-1)) # Biased estimation of our model (biais)

# R square
SCR <- sum((ldAGB-ldAGB.Mada.II.1)^2)
SCT <- sum((ldAGB-mean(ldAGB))^2)
R.square <- 1-SCR/SCT

# AIC
logL <- sum(dnorm(ldAGB,mean=ldAGB.Mada.II.1,sd=sd(Res),log=TRUE))
AIC <- -2*logL+2*npar

# Completing the result matrix
Result.Mat$beta0[Result.Mat$Model=="Mada.II.1"] <- M$coeff[1]+(sd(M$residuals)^2)/2
Result.Mat$beta1[Result.Mat$Model=="Mada.II.1"] <- M$coeff[2]
Result.Mat$beta2[Result.Mat$Model=="Mada.II.1"] <- "-"
Result.Mat$beta3[Result.Mat$Model=="Mada.II.1"] <- "-"
Result.Mat$beta4[Result.Mat$Model=="Mada.II.1"] <- M$coeff[3]
Result.Mat$df[Result.Mat$Model=="Mada.II.1"] <- df
Result.Mat$RSE[Result.Mat$Model=="Mada.II.1"] <- RSE
Result.Mat$R2[Result.Mat$Model=="Mada.II.1"] <- R.square
Result.Mat$AIC[Result.Mat$Model=="Mada.II.1"] <- AIC
Result.Mat$Bias[Result.Mat$Model=="Mada.II.1"] <- Bias

# Test
sink(file=eval(parse(text=paste("file=\"./",Levels.FType[i],"/summaryMII1",Levels.FType[i],".txt\"",sep=""))))
print(summary(M))
sink()

###############################################################
#========= Biomass-diameter regression Model II.2-WEB ========#

# Inference
Y <- ldAGB-(8/3)*lD-lrho
M <- lm(Y~1,data=Data.AGB)

# Prediction with Madagascar allometry
ldAGB.Mada.II.2.WEB <- M$coeff+(8/3)*lD+lrho

# Graph
plot(ldAGB,ldAGB.Mada.II.2.WEB,xlim=c(0,10),ylim=c(0,10))
lines(seq(0,10),seq(0,10),col="red")

# df
nobs <- length(Data.AGB$dAGB)
npar <- 1
df <- nobs-npar

# RSE
# Residuals
Res <- ldAGB-ldAGB.Mada.II.2.WEB
# Distribution of the residuals
hist(Res)
# Residuals Standard Error
RSE <- sd(Res)
# Bias
Bias.Tree.Mada.II.2.WEB <- 100*(1/exp(Res)-1)
Bias <- mean(100*(1/exp(Res)-1)) # Biased estimation of our model (biais)

# R square
SCR <- sum((ldAGB-ldAGB.Mada.II.2.WEB)^2)
SCT <- sum((ldAGB-mean(ldAGB))^2)
R.square <- 1-SCR/SCT

# AIC
logL <- sum(dnorm(ldAGB,mean=ldAGB.Mada.II.2.WEB,sd=sd(Res),log=TRUE))
AIC <- -2*logL+2*npar

# Completing the result matrix
Result.Mat$beta0[Result.Mat$Model=="Mada.II.2.WEB"] <- M$coeff[1]+(sd(M$residuals)^2)/2
Result.Mat$beta1[Result.Mat$Model=="Mada.II.2.WEB"] <- 8/3
Result.Mat$beta2[Result.Mat$Model=="Mada.II.2.WEB"] <- "-"
Result.Mat$beta3[Result.Mat$Model=="Mada.II.2.WEB"] <- "-"
Result.Mat$beta4[Result.Mat$Model=="Mada.II.2.WEB"] <- "-"
Result.Mat$df[Result.Mat$Model=="Mada.II.2.WEB"] <- df
Result.Mat$RSE[Result.Mat$Model=="Mada.II.2.WEB"] <- RSE
Result.Mat$R2[Result.Mat$Model=="Mada.II.2.WEB"] <- R.square
Result.Mat$AIC[Result.Mat$Model=="Mada.II.2.WEB"] <- AIC
Result.Mat$Bias[Result.Mat$Model=="Mada.II.2.WEB"] <- Bias

# Test
sink(file=eval(parse(text=paste("file=\"./",Levels.FType[i],"/summaryMII2",Levels.FType[i],".txt\"",sep=""))))
print(summary(M))
sink()

#######################################################################
#=====================================================================#

# Backup of the lab result matrix
# Genus
write.table(MatGenusLab,
            file=eval(parse(text=paste("file=\"./",Levels.FType[i],"/Lab-Genus-",Levels.FType[i],".txt\"",sep=""))),
            col.names=TRUE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)
# Species
write.table(MatSpeciesLab,
            file=eval(parse(text=paste("file=\"./",Levels.FType[i],"/Lab-Species-",Levels.FType[i],".txt\"",sep=""))),
            col.names=TRUE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

# Backup of the allometry result matrix
Result.Mat[Result.Mat=="-"] <- 0
Result.Mat[,c(2:6)] <- round(as.numeric(unlist(Result.Mat[,c(2:6)])),digit=3)
Result.Mat[,c(8:9)] <- round(as.numeric(unlist(Result.Mat[,c(8:9)])),digit=2)
Result.Mat[,c(10)] <- round(as.numeric(unlist(Result.Mat[,c(10)])),digit=0)
Result.Mat[,c(11)] <- round(as.numeric(unlist(Result.Mat[,c(11)])),digit=1)
Result.Mat[9,3] <- "8/3"
write.table(Result.Mat,
            file=eval(parse(text=paste("file=\"./",Levels.FType[i],"/Allometry-",Levels.FType[i],".txt\"",sep=""))),
            col.names=TRUE,
            quote=FALSE,
            sep="\t",
            row.names=FALSE)

#############################
# Graphics by forest type

if(Levels.FType[i]=="dry") {

  # Covariates
  D.seq <- seq(from=1,to=35,length.out=100)
  rho.m <- round(mean(Data.AGB$DryDensityTree,na.rm=TRUE),2)
  rho.q1 <- round(quantile(Data.AGB$DryDensityTree,0.025,na.rm=TRUE),2)
  rho.q3 <- round(quantile(Data.AGB$DryDensityTree,0.975,na.rm=TRUE),2)
  # Parameter value height-D
  a0 <- 12.120
  a1 <- 1.300
  a2 <- 0.052
  # Height predictions
  H.hat <- a0-(a0-a1)*exp(-a2*D.seq)
  # Parameter value
  b0 <- Result.Mat$beta0[Result.Mat$Model=="Mada.I.1"]
  b1 <- as.numeric(Result.Mat$beta1[Result.Mat$Model=="Mada.I.1"])
  b2 <- as.numeric(Result.Mat$beta2[Result.Mat$Model=="Mada.I.1"])
  b4 <- Result.Mat$beta4[Result.Mat$Model=="Mada.I.1"]
  # Prediction
  AGB.m <- exp(b0+b1*log(D.seq)+b2*log(H.hat)+b4*log(rho.m))
  AGB.q1 <- exp(b0+b1*log(D.seq)+b2*log(H.hat)+b4*log(rho.q1))
  AGB.q3 <- exp(b0+b1*log(D.seq)+b2*log(H.hat)+b4*log(rho.q3))

  #########################
  # Plot predicted-observed
  postscript(file="Predicted-Observed.eps",height=5,width=10,paper="special",horizontal=FALSE)
  par(cex=1.4,mar=c(4,4,0,0),mfrow=c(1,2))
  plot(Data.AGB$DBH,Data.AGB$dAGB,
       xlim=c(0,35),
       ylim=c(0,500),
       axes=FALSE,
       col=grey(0.5),
       main="",
       xlab="D (cm)",
       ylab="AGB (Mg)")
  axis(1,at=seq(0,35,5),labels=seq(0,35,5))
  axis(2,at=seq(0,500,100),labels=seq(0,0.5,0.1))
  lines(D.seq,AGB.m,lwd=2,lty=1)
  lines(D.seq,AGB.q1,lwd=1,lty=1)
  lines(D.seq,AGB.q3,lwd=3,lty=1)
  text(x=0,y=0.90*500,labels="Madagascar\nspiny dry forests",pos=4, cex=1)
  text(x=25,y=0,labels=expression(rho[1]==Q[0.025]),pos=4)
  text(x=20,y=100,labels=expression(rho[2]==textstyle(mean)),pos=4)
  text(x=12,y=250,labels=expression(rho[3]==Q[0.975]),pos=4)

  ###################
  # Plot for the bias
  b1 <- lowess(Data.AGB$DBH,Bias.Tree.Brown)$y
  b2 <- lowess(Data.AGB$DBH,Bias.Tree.Chave.H)$y
  b3 <- lowess(Data.AGB$DBH,Bias.Tree.Chave.D)$y
  b4 <- lowess(Data.AGB$DBH,Bias.Tree.Mada.I.1)$y
  b5 <- lowess(Data.AGB$DBH,Bias.Tree.Mada.I.1.wHD)$y
  b6 <- lowess(Data.AGB$DBH,Bias.Tree.Mada.II.2.WEB)$y
  min.Bias <- min(b1,b2,b3,b4,b5,b6)
  max.Bias <- max(b1,b2,b3,b4,b5,b6)

  postscript(file="Bias.eps",height=5,width=10,paper="special",horizontal=FALSE)
  par(cex=1.4,mar=c(4,4,0,0),mfrow=c(1,2))
  plot(lowess(Data.AGB$DBH,Bias.Tree.Brown),
       bty="n",
       type="l",
       col=grey(0.5),
       xlim=c(0,35),
       ylim=c(-50,200),
       xlab="D (cm)",
       ylab="Bias (% of AGB)",
       lwd=2,
       lty=5)
  segments(x0=0,y0=0,x1=35,y1=0,lwd=1,col=grey(0.7),lty=1)
  
  lines(lowess(Data.AGB$DBH,Bias.Tree.Chave.H),lty=1,lwd=2,col=grey(0.5))
  lines(lowess(Data.AGB$DBH,Bias.Tree.Chave.D),lty=4,lwd=2,col=grey(0.5))
  lines(lowess(Data.AGB$DBH,Bias.Tree.Mada.I.1),lty=1,lwd=3)
  lines(lowess(Data.AGB$DBH,Bias.Tree.Mada.I.1.wHD),lty=2,lwd=3)
  lines(lowess(Data.AGB$DBH,Bias.Tree.Mada.II.2.WEB),lty=2,lwd=2,col=grey(0.5))
  
  text(x=0,y=180,labels="Madagascar\nspiny dry forests",pos=4,cex=1)
  
  legend(x=0,y=170,legend=Result.Mat$Model[c(1:5,9)],
         lty=c(5,1,4,1,2,2),
         lwd=c(2,2,2,3,3,2),
         col=c(rep(grey(0.5),3),rep("black",2),grey(0.5)),
         bty="n",cex=0.8)

  ##########################
  # Plot for the predictions

  # rho.m
  #rho.m <- 0.10
  
  # Max DBH and x.seq
  x.seq.GOP <- seq(from=0,to=48,length.out=100)
  x.seq.Brown <- seq(from=0,to=30,length.out=100)
  x.seq.Chave <- seq(from=0,to=63.4,length.out=100)
  x.seq.WEB <- seq(from=0,to=100,length.out=100)

  # Parameter value
  b0 <- Result.Mat$beta0[Result.Mat$Model=="Mada.I.1"]
  b1 <- as.numeric(Result.Mat$beta1[Result.Mat$Model=="Mada.I.1"])
  b2 <- as.numeric(Result.Mat$beta2[Result.Mat$Model=="Mada.I.1"])
  b4 <- Result.Mat$beta4[Result.Mat$Model=="Mada.I.1"]

  # Predictions
  H.seq.GOP <- 12.12-(12.12-1.3)*exp(-0.052*x.seq.GOP)
  y.seq.GOP <- exp(b0+b1*log(x.seq.GOP)+b2*log(H.seq.GOP)+b4*log(rho.m))
  y.seq.Brown <- exp(-1.473+2*log(x.seq.Brown))
  y.seq.Chave <- exp(-0.667+1.784*log(x.seq.Chave)+
                     0.207*log(x.seq.Chave)^2-0.028*log(x.seq.Chave)^3+log(rho.m))
  y.seq.WEB <- exp(-1.999+(8/3)*log(x.seq.WEB)+log(rho.m))

  # Plot
  
  postscript(file="Pred-BigTrees.eps",height=5,width=10,paper="special",horizontal=FALSE)
  par(cex=1.4,mar=c(4,4,0,0),mfrow=c(1,2))
  plot(x.seq.WEB,y.seq.WEB,
       xlim=c(0,70),
       ylim=c(0,3000),
       xlab="D (cm)",
       ylab="AGB (Mg)",
       axes=FALSE,
       lty=2,lwd=2,col=grey(0.5),
       type="l")
  axis(1,at=seq(0,70,10),labels=seq(0,70,10))
  axis(2,at=seq(0,3000,500),labels=seq(0,3,0.5))
  lines(x.seq.Brown,y.seq.Brown,lty=5,lwd=2,col=grey(0.5))
  lines(x.seq.Chave,y.seq.Chave,lty=4,lwd=2,col=grey(0.5))
  lines(x.seq.GOP,y.seq.GOP,lty=2,lwd=3)
  legend(x=0,y=0.80*3000,legend=Result.Mat$Model[c(1,3,5,9)],
         lty=c(5,4,2,2),
         lwd=c(2,2,3,2),
         col=c(rep(grey(0.5),2),"black",grey(0.5)),
         bty="n",cex=0.8)
  text(x=0,y=0.5*3000,labels=expression(rho==0.56),cex=0.9,pos=4)
  text(x=0,y=0.90*3000,labels="Madagascar\nspiny dry forests",pos=4,cex=1)

  dev.new()
}



if(Levels.FType[i]=="moist-wet") {

  # Covariates
  D.seq <- seq(from=1,to=60,length.out=100)
  rho.m <- round(mean(Data.AGB$DryDensityTree,na.rm=TRUE),2)
  rho.q1 <- round(quantile(Data.AGB$DryDensityTree,0.025,na.rm=TRUE),2)
  rho.q3 <- round(quantile(Data.AGB$DryDensityTree,0.975,na.rm=TRUE),2)
  # Parameter value
  b0 <- Result.Mat$beta0[Result.Mat$Model=="Mada.II.1"]
  b1 <- as.numeric(Result.Mat$beta1[Result.Mat$Model=="Mada.II.1"])
  b4 <- Result.Mat$beta4[Result.Mat$Model=="Mada.II.1"]
  # Prediction
  AGB.m <- exp(b0+b1*log(D.seq)+b4*log(rho.m))
  AGB.q1 <- exp(b0+b1*log(D.seq)+b4*log(rho.q1))
  AGB.q3 <- exp(b0+b1*log(D.seq)+b4*log(rho.q3))

  #########################
  # Plot predicted-observed
  dev.set(3)
  plot(Data.AGB$DBH,Data.AGB$dAGB,
       xlim=c(0,60),
       ylim=c(0,2500),
       axes=FALSE,
       col=grey(0.5),
       main="",
       xlab="D (cm)",
       ylab="AGB (Mg)")
  axis(1,at=seq(0,60,10),labels=seq(0,60,10))
  axis(2,at=seq(0,2500,500),labels=seq(0,2.5,0.5))
  lines(D.seq,AGB.m,lwd=2,lty=1)
  lines(D.seq,AGB.q1,lwd=1,lty=1)
  lines(D.seq,AGB.q3,lwd=3,lty=1)
  text(x=0,y=0.90*2500,labels=paste("Madagascar\n",Levels.FType[i]," forests",sep=""),pos=4, cex=1)
  text(x=40,y=600,labels=expression(rho[1]==Q[0.025]),pos=4)
  text(x=24,y=1450,labels=expression(rho[2]==textstyle(mean)),pos=4)
  text(x=26,y=1650,labels=expression(rho[3]==Q[0.975]),pos=4)
  arrows(x0=35,y0=1400,x1=D.seq[69],y1=AGB.m[69],length=0.1)

  ###################
  # Plot for the bias
  b1 <- lowess(Data.AGB$DBH,Bias.Tree.Brown)$y
  b2 <- lowess(Data.AGB$DBH,Bias.Tree.Chave.H)$y
  b3 <- lowess(Data.AGB$DBH,Bias.Tree.Chave.D)$y
  b4 <- lowess(Data.AGB$DBH,Bias.Tree.Mada.I.1)$y
  b5 <- lowess(Data.AGB$DBH,Bias.Tree.Mada.II.1)$y
  b6 <- lowess(Data.AGB$DBH,Bias.Tree.Mada.II.2.WEB)$y
  min.Bias <- min(b1,b2,b3,b4,b5,b6)
  max.Bias <- max(b1,b2,b3,b4,b5,b6)
  
  dev.set(4)
  plot(lowess(Data.AGB$DBH,Bias.Tree.Brown),
       bty="n",
       axes=FALSE,
       col=grey(0.5),
       type="l",
       xlim=c(0,60),
       ylim=c(-40,60),
       xlab="D (cm)",
       ylab="Bias (% of AGB)",
       lwd=2,
       lty=5)
  axis(1,at=seq(from=0,to=60,by=10),labels=seq(from=0,to=60,by=10))
  axis(2,at=seq(from=-40,to=60,by=20),labels=seq(from=-40,to=60,by=20))
  segments(x0=0,y0=0,x1=60,y1=0,lwd=1,col=grey(0.7),lty=1)

  lines(lowess(Data.AGB$DBH,Bias.Tree.Chave.H),lty=1,lwd=2,col=grey(0.5))
  lines(lowess(Data.AGB$DBH,Bias.Tree.Chave.D),lty=4,lwd=2,col=grey(0.5))
  lines(lowess(Data.AGB$DBH,Bias.Tree.Mada.I.1),lty=1,lwd=3)
  lines(lowess(Data.AGB$DBH,Bias.Tree.Mada.II.1),lty=2,lwd=3)
  lines(lowess(Data.AGB$DBH,Bias.Tree.Mada.II.2.WEB),lty=2,lwd=2,col=grey(0.5))
  
  text(x=0,y=52,labels=paste("Madagascar\n",Levels.FType[i]," forests",sep=""),pos=4,cex=1)

  legend(x=20,y=-10,legend=Result.Mat$Model[c(1:4,8,9)],
         lty=c(5,1,4,1,2,2),
         lwd=c(2,2,2,3,3,2),
         col=c(rep(grey(0.5),3),rep("black",2),grey(0.5)),
         bty="n",cex=0.8)

  ##########################
  # Plot for the predictions

  # Max DBH and x.seq
  x.seq.GOP <- seq(from=0,to=130,length.out=100)
  x.seq.Brown <- seq(from=0,to=148,length.out=100)
  x.seq.Chave <- seq(from=0,to=156,length.out=100)
  x.seq.WEB <- seq(from=0,to=200,length.out=100)

  # Parameter value
  b0 <- Result.Mat$beta0[Result.Mat$Model=="Mada.II.1"]
  b1 <- as.numeric(Result.Mat$beta1[Result.Mat$Model=="Mada.II.1"])
  b4 <- Result.Mat$beta4[Result.Mat$Model=="Mada.II.1"]

  # Predictions
  y.seq.GOP <- exp(b0+b1*log(x.seq.GOP)+b4*log(rho.m))
  y.seq.Brown <- exp(-2.289+2.649*log(x.seq.Brown)-0.021*log(x.seq.Brown)^2)
  y.seq.Chave <- exp(-1.499+2.148*log(x.seq.Chave)+
                     0.207*log(x.seq.Chave)^2-0.028*log(x.seq.Chave)^3+log(rho.m))
  y.seq.WEB <- exp(-2.122+(8/3)*log(x.seq.WEB)+log(rho.m))

  # Plot
  dev.set(5)
  plot(x.seq.WEB,y.seq.WEB,
       xlim=c(0,200),
       ylim=c(0,50000),
       xlab="D (cm)",
       ylab="AGB (Mg)",
       axes=FALSE,
       lty=2,lwd=2,col=grey(0.5),
       type="l")
  axis(1,at=seq(0,200,40),labels=seq(0,200,40))
  axis(2,at=seq(0,50000,10000),labels=seq(0,50,10))
  lines(x.seq.Brown,y.seq.Brown,lty=5,lwd=2,col=grey(0.5))
  lines(x.seq.Chave,y.seq.Chave,lty=4,lwd=2,col=grey(0.5))
  lines(x.seq.GOP,y.seq.GOP,lty=2,lwd=3)
  legend(x=0,y=0.80*50000,legend=Result.Mat$Model[c(1,3,8,9)],
         lty=c(5,4,2,2),
         lwd=c(2,2,3,2),
         col=c(rep(grey(0.5),2),"black",grey(0.5)),
         bty="n",cex=0.8)
  text(x=0,y=0.90*50000,labels=paste("Madagascar\n",Levels.FType[i]," forests",sep=""),pos=4,cex=1)
  text(x=0,y=0.5*50000,labels=expression(rho==0.59),cex=0.9,pos=4)

  graphics.off()
}

}


##########################################################################################################################
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
##
#                                            End
##
#/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
##########################################################################################################################
