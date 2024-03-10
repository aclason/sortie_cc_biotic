####This script contains the code for parameterizing SORTIE for climate and nitrogen sensitive growth
#code was started by Will and Kiri Daust, July 2018
###adjusted by Erica Lilles October 2018 to fit SIBEC data
###redone by Erica Feb 2024 to parameterize SORTIE using SIBEC data converted to "MaxGrowth"
rm(list=ls())

###Model aSMR <-> rSMR
library(tcltk)
library(reshape2)
library(foreach)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(data.table)
library(gridExtra)
library(mgcv)
library(rgl)
library(deldir)
library(likelihood)
library(lattice)

###Create rSMR -> aSMR crosswalk######################
setwd("C:/Users/elilles/OneDrive - Government of BC/SIBEC likelihood")

####Create SI models for each species####
###__________________________________________________________________#####
load ("rSMR_aSMR_CalcList.RData")
Eda <- read.csv("Edatopic_v10.7.csv")

colnames(SMRCross) <- c("BGC", "rSMR", "aSMR")
SMRCross$rSMR <- gsub("[[:alpha:]]","", SMRCross$rSMR)
SMRCross$aSMR <- round(SMRCross$aSMR, digits = 0)

sibec <- read.csv("SIBEC_for_Portfolio2.csv") ###import actual SI values
sibec<-subset(sibec, sibec$PlotCountSpp>0 )## remove site series with no plots

#which species have enough data to make predictive model? only Sx, Fd, Pl and Bl have over 100
length(sibec$TreeSpp[sibec$TreeSpp == "Sx"]) #215
length(sibec$TreeSpp[sibec$TreeSpp == "Fd"]) #134
length(sibec$TreeSpp[sibec$TreeSpp == "Pl"]) #236
length(sibec$TreeSpp[sibec$TreeSpp == "Sw"]) #0
length(sibec$TreeSpp[sibec$TreeSpp == "Se"]) #28
length(sibec$TreeSpp[sibec$TreeSpp == "Ss"]) #22
length(sibec$TreeSpp[sibec$TreeSpp == "Cw"]) #46
length(sibec$TreeSpp[sibec$TreeSpp == "Lw"]) #14
length(sibec$TreeSpp[sibec$TreeSpp == "Bl"]) #125
length(sibec$TreeSpp[sibec$TreeSpp == "Ba"]) #22
length(sibec$TreeSpp[sibec$TreeSpp == "Hw"]) #58
length(sibec$TreeSpp[sibec$TreeSpp == "Hm"]) #0
length(sibec$TreeSpp[sibec$TreeSpp == "Yc"]) #0

##############below code has already been used to prep data for each species. It needs to be
#re-examined in future iterations
sibec <- sibec[sibec$TreeSpp == "Pl",]##choose species
sibec <- sibec[!is.na(sibec$BGCUnit),]
#Erica note: as is, this is giving the number of site series listed in "sibec" for each BGCUnit
#I don't think this is what we want because MeanPlotSiteIndex also has si values with no plots listed, which must be guesses?
#numPlots <- aggregate(MeanPlotSiteIndex ~ BGCUnit, sibec, FUN = length)
#BGC <- numPlots$BGCUnit[numPlots$MeanPlotSiteIndex > 4]###remove SI values in zones with little data (use 2 for Lw and Py)
#sibec <- sibec[sibec$BGCUnit %in% BGC,]
#instead, lets first remove site series that don't have any plots...this is done above in the code
sibec_plots<-subset(sibec, sibec$PlotCountSpp>0 )
min(sibec_plots$PlotCountSpp) #remaining site series have at least 7 plots with data (for Bl), that seems sufficient to include
#then we can use above code but change numPlots to numSiteSeries which describes the object more accurately
numSiteSeries <- aggregate(SiteSeries ~ BGCUnit, sibec_plots, FUN = length)
BGC <- numSiteSeries$BGCUnit[numSiteSeries$SiteSeries >= 4] #changed to >= 4 for Bl so enough BGCs would be included
sibec <- sibec_plots[sibec_plots$BGCUnit %in% BGC,]

####for BL, 3 BGC  unit site series don't match up between sibec data and edatope so we are losing data during this merge
#### need another cross walk table here for "ESSFdc1" "ESSFwc1" and "SBSwk3"(doesn't have data for edatopic for SBSwk3/08 and another site series)
#SBSwk3 doesn't get through the merge while looped but it does individually because the Edatopic is missing for 2 site series and there is a is.na function

#try converting sibec$BGCUnit to character to fix loop error
sibec$BGCUnit<-as.character(sibec$BGCUnit)
Eda$MergedBGC<-as.character(Eda$MergedBGC)
###Loop to match SI values with edatopes and convert to aSMR
out <- foreach(i = unique(sibec$BGCUnit), .combine = rbind) %do% {
  SIsub <- sibec[sibec$BGCUnit == i,c(9,7)] #extract SS_nospace and mean site index
  edaSub <- Eda[Eda$MergedBGC == i, 3:4] #extract SS_nospace and edatopic position
  if(length(edaSub$SS_NoSpace) > 0){
    SIsub <- merge(SIsub, edaSub, by = "SS_NoSpace", all.x = TRUE) #merge site index with edatopic position
    if(!any(is.na(SIsub$Edatopic))){
      SIeda <- aggregate(MeanPlotSiteIndex ~ Edatopic, data =  SIsub, FUN = mean) #this combines duplicates of site series and edatopic position when they are listed in more than one region
      colnames(SIeda)[2] <- "SI"
      SIeda$rSMR <- gsub("[[:alpha:]]","", SIeda$Edatopic) #separate SMR from Edatopic position
      SIeda$SNR <- gsub("[[:digit:]]","", SIeda$Edatopic)  #separate SNR from Edatopic position
      SIeda$rSMR <- as.numeric(as.character(SIeda$rSMR))
      SMRsub <- SMRCross[SMRCross$BGC == i,-1] #get aSMR rSMR conversion table
      SIeda <- merge(SIeda, SMRsub, by = "rSMR", all.x = TRUE) #merge in aSMR
      SIeda<-unique(SIeda) #remove duplicated rows that occured when a aSMR class included more than one rSMR
      #Erica note: taking out below commands because I prefer to work with the raw values (not standardized), at least until I understand better
      #temp <- SIeda[SIeda$aSMR %% 1 != 0,] ###if 0.5, give to higher and lower whole num
      #temp$aSMR <- ceiling(temp$aSMR)
      #SIeda$aSMR <- floor(SIeda$aSMR)
      #SIeda <- rbind(SIeda,temp)

      #SIeda$SI <- SIeda$SI/max(SIeda$SI)###Standardise out of max SI value in BGC
      SIeda <- SIeda[,c(4,5,3)]
      SIeda$BGC <- i
      SIeda
    }}}
SI2 <- out

#New exploratory code by Erica
head(SI2)
SI2$SNR<-ifelse(SI2$SNR == "A", 1, ifelse(SI2$SNR == "B", 2, ifelse(SI2$SNR == "C", 3,ifelse(SI2$SNR == "D", 4, 5))))
#write data to file so above steps can be skipped in future
write.csv(SI2, "SiteIndex_forLikelihood_Pl.csv")



#########################################################################################
#load data that was prepped by above code for each species
SI2_Pl<-read.csv("SiteIndex_forLikelihood_Pl.csv")
SI2_Bl<-read.csv("SiteIndex_forLikelihood_Bl.csv")
SI2_Sx<-read.csv("SiteIndex_forLikelihood_Sx.csv")
SI2_At<-read.csv("SiteIndex_forLikelihood_At.csv")

################################
################################ combine with climate data and try model with climate vars affecting maxSI
################################

dat <- fread("BECv11_100Pt_Normal_1961_1990MSY.csv", data.table = FALSE)
dat <- dat[,c("ID2","Elevation", "DD5", "CMD" )]
dat$ID2 <- gsub("[[:space:]]","",dat$ID2)

mean_dat<-ddply(dat, .(ID2), numcolwise(mean)) #average the 100 points for each BGC for elevation and all climate variables

SI_climate_Pl<-merge(SI2_Pl, mean_dat, by.x = "BGC", by.y = "ID2", all.x = TRUE)

plot(SI_climate_Pl$DD5, SI_climate_Pl$SI)
plot(SI_climate_Pl$CMD, SI_climate_Pl$SI)

SI_climate_Bl<-merge(SI2_Bl, mean_dat, by.x = "BGC", by.y = "ID2", all.x = TRUE)
SI_climate_Sx<-merge(SI2_Sx, mean_dat, by.x = "BGC", by.y = "ID2", all.x = TRUE)
SI_climate_At<-merge(SI2_At, mean_dat, by.x = "BGC", by.y = "ID2", all.x = TRUE)

###################################################################################################################
###################################################################################################################
####################################         Pine            ######################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

#######################
############################# SORTIE aSMR double logistic; double logistic for temp and including SNR
###################### fitting Pine without ICHwk4 or MSdm1 which give it a big jog... sleuthing required to figure out if climate data is wrong or if a variable is missing
###################### warning this parameterization goes steeply down at high DD5s so is not suitable for hotter climates until there is more data on the high end
SI_climate_Pl.noICHwk4<-subset(SI_climate_Pl, SI_climate_Pl$BGC !=  "ICHwk4")
SI_climate_Pl.noICHwk4.noMSdm1<-subset(SI_climate_Pl.noICHwk4, SI_climate_Pl.noICHwk4$BGC !=  "MSdm1")

SORTIE.aSMRSNR.temp.model.log2 <- function(maxSI, aL, bL, cL, aH, bH, cH, aSMR, dL, eL, fL, dH, eH, fH,  temp, Xo, Xb, SNR)
{  maxSI * (aL+ (1-aL)/(1+ ((bL/aSMR)^cL))) * (aH+ (1-aH)/(1+ ((aSMR/bH)^cH)) ) * (dL+ (1-dL)/(1+ ((eL/temp)^fL))) * (dH+ (1-dH)/(1+ ((temp/eH)^fH)) ) * exp(-0.5 * ((SNR - Xo)/Xb)^2) }

var <- list(aSMR = "aSMR", temp = "DD5", SNR = "SNR" )
par <- list( maxSI = 42, aL = 0.001, bL = 3, cL = .5, aH =.05, bH = 8, cH = 70, dL = 0.1, eL = 500, fL = 1, dH =4, eH = 500, fH = 0.01,  Xo = 3, Xb = 1, sd = 1)
par_lo <- list( maxSI = 0,  aL = 0,     bL = 0,  cL = 0, aH =0, bH = 0,  cH = 0, dL = 0,  eL = 0,   fL = 0,  dH =0,  eH = 0,    fH = 0, Xo = 0, Xb = 0,    sd = 0)
par_hi <- list( maxSI = 150, aL = .01, bL = 6, cL = 10, aH =.1, bH = 15, cH = 100, dL = 1, eL = 2000, fL = 20, dH =1, eH = 2000, fH = 500, Xo = 50, Xb = 35,   sd = 10)

var$x <- "SI" #this is the site index
var$mean <- "predicted"
var$log <- TRUE

model_results_Pl <- anneal (model = SORTIE.aSMRSNR.temp.model.log2, par = par, var = var, source_data = SI_climate_Pl.noICHwk4.noMSdm1, par_lo = par_lo, par_hi = par_hi,
                               pdf=dnorm, dep_var="SI",   initial_temp = 5, temp_red = 0.96, max_iter=100000, hessian = FALSE)

write_results(model_results_Pl,"SORTIE parameterization/Pl1_SiteIndex.txt", data=F)

#look at residuals
plot(model_results_Pl$source_data$predicted, model_results_Pl$source_data$SI)
SBSmc2_points<-subset(model_results_Pl$source_data, model_results_Pl$source_data$BGC == "SBSmc2")
points(SBSmc2_points$predicted, SBSmc2_points$SI, pch = 2, col = "red")
abline(0,1)
xyplot(predicted~ SI|BGC, data = model_results_Pl$source_data)
xyplot(predicted~ SI, group = BGC, data = model_results_Pl$source_data)

resid<- model_results_Pl$source_data$predicted-model_results_Pl$source_data$SI
plot(model_results_Pl$source_data$SI, resid)
##  plot a 0 line  (i.e. intercept = 0, slope = 0)
abline(0,0)

plot(model_results_Pl$source_data$aSMR, model_results_Pl$source_data$SI)
x<-seq(1, 8, .1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Pl$best_pars$maxSI , aL = model_results_Pl$best_pars$aL, bL = model_results_Pl$best_pars$bL, cL = model_results_Pl$best_pars$cL,
                               aH = model_results_Pl$best_pars$aH, bH = model_results_Pl$best_pars$bH, cH = model_results_Pl$best_pars$cH,
                               dL = model_results_Pl$best_pars$dL, eL = model_results_Pl$best_pars$eL, fL = model_results_Pl$best_pars$fL,
                               dH = model_results_Pl$best_pars$dH, eH = model_results_Pl$best_pars$eH, fH = model_results_Pl$best_pars$fH,
                               Xo = model_results_Pl$best_pars$Xo, Xb = model_results_Pl$best_pars$Xb,
                               SNR = mean(model_results_Pl$source_data$SNR), aSMR = x, temp = mean(model_results_Pl$source_data$DD5))
lines(x,y)

plot(model_results_Pl$source_data$DD5, model_results_Pl$source_data$SI, xlim = c(300, 1500))
x<-seq(300, 1500, 1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Pl$best_pars$maxSI , aL = model_results_Pl$best_pars$aL, bL = model_results_Pl$best_pars$bL, cL = model_results_Pl$best_pars$cL,
                               aH = model_results_Pl$best_pars$aH, bH = model_results_Pl$best_pars$bH, cH = model_results_Pl$best_pars$cH,
                               dL = model_results_Pl$best_pars$dL, eL = model_results_Pl$best_pars$eL, fL = model_results_Pl$best_pars$fL,
                               dH = model_results_Pl$best_pars$dH, eH = model_results_Pl$best_pars$eH, fH = model_results_Pl$best_pars$fH,
                               Xo = model_results_Pl$best_pars$Xo, Xb = model_results_Pl$best_pars$Xb,
                               SNR = mean(model_results_Pl$source_data$SNR), temp = x, aSMR = mean(model_results_Pl$source_data$aSMR))
lines(x,y)
SBSmc2_points<-subset(SI_climate_Pl, SI_climate_Pl$BGC == "SBSmc2")
points(SBSmc2_points$DD5, SBSmc2_points$SI, pch = 2, col = "red")
SI_climate_Pl[which(SI_climate_Pl$DD5 >1050),]

plot(model_results_Pl$source_data$SNR, model_results_Pl$source_data$SI)
x<-seq(1, 5, .1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Pl$best_pars$maxSI , aL = model_results_Pl$best_pars$aL, bL = model_results_Pl$best_pars$bL, cL = model_results_Pl$best_pars$cL,
                                  aH = model_results_Pl$best_pars$aH, bH = model_results_Pl$best_pars$bH, cH = model_results_Pl$best_pars$cH,
                                  dL = model_results_Pl$best_pars$dL, eL = model_results_Pl$best_pars$eL, fL = model_results_Pl$best_pars$fL,
                                  dH = model_results_Pl$best_pars$dH, eH = model_results_Pl$best_pars$eH, fH = model_results_Pl$best_pars$fH,
                                  Xo = model_results_Pl$best_pars$Xo, Xb = model_results_Pl$best_pars$Xb,
                                  SNR = x, temp = mean(model_results_Pl$source_data$DD5), aSMR = mean(model_results_Pl$source_data$aSMR))
lines(x,y)
SBSmc2_points<-subset(SI_climate_Pl, SI_climate_Pl$BGC == "SBSmc2")
points(SBSmc2_points$SNR, SBSmc2_points$SI, pch = 2, col = "red")

###################
####################
####################
#take above function and predict new MaxGrowth parameters
###################
####################
####################
subset(SMRCross, SMRCross$BGC == "SBSmc2")

#Create conversion from Site Index to MaxGrowth
SI_conversion<-matrix(nrow = 4,ncol = 2)

#predict site index for SBSmc2 '09' sites...
SI_conversion[4,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Pl$best_pars$maxSI , aL = model_results_Pl$best_pars$aL, bL = model_results_Pl$best_pars$bL, cL = model_results_Pl$best_pars$cL,
                               aH = model_results_Pl$best_pars$aH, bH = model_results_Pl$best_pars$bH, cH = model_results_Pl$best_pars$cH,
                               dL = model_results_Pl$best_pars$dL, eL = model_results_Pl$best_pars$eL, fL = model_results_Pl$best_pars$fL,
                               dH = model_results_Pl$best_pars$dH, eH = model_results_Pl$best_pars$eH, fH = model_results_Pl$best_pars$fH,
                               Xo = model_results_Pl$best_pars$Xo, Xb = model_results_Pl$best_pars$Xb,
                               SNR = 5, temp = 875.73, aSMR = 6)
#predict site index for SBSmc2 '06' sites
SI_conversion[3,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Pl$best_pars$maxSI , aL = model_results_Pl$best_pars$aL, bL = model_results_Pl$best_pars$bL, cL = model_results_Pl$best_pars$cL,
                               aH = model_results_Pl$best_pars$aH, bH = model_results_Pl$best_pars$bH, cH = model_results_Pl$best_pars$cH,
                               dL = model_results_Pl$best_pars$dL, eL = model_results_Pl$best_pars$eL, fL = model_results_Pl$best_pars$fL,
                               dH = model_results_Pl$best_pars$dH, eH = model_results_Pl$best_pars$eH, fH = model_results_Pl$best_pars$fH,
                               Xo = model_results_Pl$best_pars$Xo, Xb = model_results_Pl$best_pars$Xb,
                               SNR = 4, temp = 875.73, aSMR = 5)
#predict site index for SBSmc2 '01' sites
SI_conversion[2,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Pl$best_pars$maxSI , aL = model_results_Pl$best_pars$aL, bL = model_results_Pl$best_pars$bL, cL = model_results_Pl$best_pars$cL,
                               aH = model_results_Pl$best_pars$aH, bH = model_results_Pl$best_pars$bH, cH = model_results_Pl$best_pars$cH,
                               dL = model_results_Pl$best_pars$dL, eL = model_results_Pl$best_pars$eL, fL = model_results_Pl$best_pars$fL,
                               dH = model_results_Pl$best_pars$dH, eH = model_results_Pl$best_pars$eH, fH = model_results_Pl$best_pars$fH,
                               Xo = model_results_Pl$best_pars$Xo, Xb = model_results_Pl$best_pars$Xb,
                               SNR = 3, temp = 875.73, aSMR = 4)

#predict site index for SBSmc2 '02' sites
SI_conversion[1,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Pl$best_pars$maxSI , aL = model_results_Pl$best_pars$aL, bL = model_results_Pl$best_pars$bL, cL = model_results_Pl$best_pars$cL,
                               aH = model_results_Pl$best_pars$aH, bH = model_results_Pl$best_pars$bH, cH = model_results_Pl$best_pars$cH,
                               dL = model_results_Pl$best_pars$dL, eL = model_results_Pl$best_pars$eL, fL = model_results_Pl$best_pars$fL,
                               dH = model_results_Pl$best_pars$dH, eH = model_results_Pl$best_pars$eH, fH = model_results_Pl$best_pars$fH,
                               Xo = model_results_Pl$best_pars$Xo, Xb = model_results_Pl$best_pars$Xb,
                               SNR = 2, temp = 875.73, aSMR = 3)

#Adult MaxGrowth for lodgepole pine is 5.4 (mm/yr) * (SNR/6)^0.64
SI_conversion[1,2]<- 5.4 * (2/5)^0.64
SI_conversion[2,2]<- 5.4 * (3/5)^0.64
SI_conversion[3,2]<- 5.4 * (4/5)^0.64
SI_conversion[4,2]<- 5.4 * (5/5)^0.64
colnames(SI_conversion)<- c("SI", "MaxGrowth")
SI_conversion<-as.data.frame(SI_conversion)
plot(SI_conversion$SI, SI_conversion$MaxGrowth)
slope<- (SI_conversion$MaxGrowth[4]-SI_conversion$MaxGrowth[1])/(SI_conversion$SI[4]-SI_conversion$SI[1])
# y = mx + b
SI_conversion$MaxGrowth[4] = SI_conversion$SI[4] * slope + b
SI_conversion$MaxGrowth[4]-SI_conversion$SI[4] * slope
slope
SI_conv_fnct<-function(x) {0.6217102 * x + -7.384596}
x<- seq(13, 21, .1)
y<- SI_conv_fnct(x)
lines(x, y)

SI_climate_Pl.noICHwk4.noMSdm1$MaxGrowth_conv<- SI_conv_fnct(SI_climate_Pl.noICHwk4.noMSdm1$SI)

##################                                ####################
################## relative productivity fitting ###################
#################                                 ###################

SORTIE.climate.sens <- function(MaxGrowth, aL, bL, cL, aH, bH, cH, aSMR, dL, eL, fL, dH, eH, fH,  temp, Xo, Xb, SNR)
{  MaxGrowth * (aL+ (1-aL)/(1+ ((bL/aSMR)^cL))) * (aH+ (1-aH)/(1+ ((aSMR/bH)^cH)) ) * (dL+ (1-dL)/(1+ ((eL/temp)^fL))) * (dH+ (1-dH)/(1+ ((temp/eH)^fH)) ) * exp(-0.5 * ((SNR - Xo)/Xb)^2) }

var <- list(aSMR = "aSMR", temp = "DD5", SNR = "SNR" )
par <- list( MaxGrowth = 1, aL = 0.001, bL = 3, cL = .5, aH =.05, bH = 8, cH = 70, dL = 0.1, eL = 500, fL = 1, dH =4, eH = 500, fH = 0.01,  Xo = 3, Xb = 1, sd = 1)
par_lo <- list( MaxGrowth = 0,  aL = 0,     bL = 0,  cL = 0, aH =0, bH = 0,  cH = 0, dL = 0,  eL = 0,   fL = 0,  dH =0,  eH = 0,    fH = 0, Xo = 0, Xb = 0,    sd = 0)
par_hi <- list( MaxGrowth = 60, aL = .01, bL = 6, cL = 10, aH =.1, bH = 15, cH = 50, dL = 1, eL = 2000, fL = 20, dH =1, eH = 2000, fH = 500, Xo = 40, Xb = 35,   sd = 10)

var$x <- "MaxGrowth_conv" #this is the site index
var$mean <- "predicted"
var$log <- TRUE

model5_results_Pl <- anneal (model = SORTIE.climate.sens, par = par, var = var, source_data = SI_climate_Pl.noICHwk4.noMSdm1, par_lo = par_lo, par_hi = par_hi,
                               pdf=dnorm, dep_var="MaxGrowth_conv",   initial_temp = 5, temp_red = 0.96, max_iter=100000, hessian = FALSE)

write_results(model5_results_Pl,"SORTIE parameterization/Pl3_MaxGrowth_adult.txt", data=F)

#look at residuals
plot(model5_results_Pl$source_data$predicted, model5_results_Pl$source_data$MaxGrowth_conv)
SBSmc2_points<-subset(model5_results_Pl$source_data, model5_results_Pl$source_data$BGC == "SBSmc2")
points(SBSmc2_points$predicted, SBSmc2_points$MaxGrowth_conv, pch = 2, col = "red")
abline(0,1)
xyplot(predicted~ MaxGrowth_conv|BGC, data = model5_results_Pl$source_data)
xyplot(predicted~ MaxGrowth_conv, group = BGC, data = model5_results_Pl$source_data)

resid<- model5_results_Pl$source_data$predicted-model5_results_Pl$source_data$MaxGrowth_conv
plot(model5_results_Pl$source_data$MaxGrowth_conv, resid)
##  plot a 0 line  (i.e. intercept = 0, slope = 0)
abline(0,0)

plot(model5_results_Pl$source_data$aSMR, model5_results_Pl$source_data$MaxGrowth_conv)
x<-seq(1, 8, .1)
y<-SORTIE.climate.sens(MaxGrowth=model5_results_Pl$best_pars$MaxGrowth , aL = model5_results_Pl$best_pars$aL, bL = model5_results_Pl$best_pars$bL, cL = model5_results_Pl$best_pars$cL,
                                  aH = model5_results_Pl$best_pars$aH, bH = model5_results_Pl$best_pars$bH, cH = model5_results_Pl$best_pars$cH,
                                  dL = model5_results_Pl$best_pars$dL, eL = model5_results_Pl$best_pars$eL, fL = model5_results_Pl$best_pars$fL,
                                  dH = model5_results_Pl$best_pars$dH, eH = model5_results_Pl$best_pars$eH, fH = model5_results_Pl$best_pars$fH,
                                  Xo = model5_results_Pl$best_pars$Xo, Xb = model5_results_Pl$best_pars$Xb,
                                  SNR = mean(model5_results_Pl$source_data$SNR), aSMR = x, temp = mean(model5_results_Pl$source_data$DD5))
lines(x,y)

plot(model5_results_Pl$source_data$DD5, model5_results_Pl$source_data$MaxGrowth_conv)
x<-seq(300, 1500, 1)
y<-SORTIE.climate.sens(MaxGrowth=model5_results_Pl$best_pars$MaxGrowth , aL = model5_results_Pl$best_pars$aL, bL = model5_results_Pl$best_pars$bL, cL = model5_results_Pl$best_pars$cL,
                                  aH = model5_results_Pl$best_pars$aH, bH = model5_results_Pl$best_pars$bH, cH = model5_results_Pl$best_pars$cH,
                                  dL = model5_results_Pl$best_pars$dL, eL = model5_results_Pl$best_pars$eL, fL = model5_results_Pl$best_pars$fL,
                                  dH = model5_results_Pl$best_pars$dH, eH = model5_results_Pl$best_pars$eH, fH = model5_results_Pl$best_pars$fH,
                                  Xo = model5_results_Pl$best_pars$Xo, Xb = model5_results_Pl$best_pars$Xb,
                                  SNR = mean(model5_results_Pl$source_data$SNR), temp = x, aSMR = mean(model5_results_Pl$source_data$aSMR))
lines(x,y)
#can a new maxgrowth be put in site index function? Answer is No
SI_conv_fnct(model5_results_Pl$best_pars$MaxGrowth)
y2<-SORTIE.aSMRSNR.temp.model.log2(maxSI=23.7  , aL = model_results_Pl$best_pars$aL, bL = model_results_Pl$best_pars$bL, cL = model_results_Pl$best_pars$cL,
                                  aH = model_results_Pl$best_pars$aH, bH = model_results_Pl$best_pars$bH, cH = model_results_Pl$best_pars$cH,
                                  dL = model_results_Pl$best_pars$dL, eL = model_results_Pl$best_pars$eL, fL = model_results_Pl$best_pars$fL,
                                  dH = model_results_Pl$best_pars$dH, eH = model_results_Pl$best_pars$eH, fH = model_results_Pl$best_pars$fH,
                                  Xo = model_results_Pl$best_pars$Xo, Xb = model_results_Pl$best_pars$Xb,
                                  SNR = mean(model5_results_Pl$source_data$SNR), temp = x, aSMR = mean(model5_results_Pl$source_data$aSMR))
lines(x,y2, lwd = 2, col = "red")


plot(model5_results_Pl$source_data$SNR, model5_results_Pl$source_data$MaxGrowth_conv)
x<-seq(1, 5, .1)
y<-SORTIE.climate.sens(MaxGrowth=model5_results_Pl$best_pars$MaxGrowth , aL = model5_results_Pl$best_pars$aL, bL = model5_results_Pl$best_pars$bL, cL = model5_results_Pl$best_pars$cL,
                                  aH = model5_results_Pl$best_pars$aH, bH = model5_results_Pl$best_pars$bH, cH = model5_results_Pl$best_pars$cH,
                                  dL = model5_results_Pl$best_pars$dL, eL = model5_results_Pl$best_pars$eL, fL = model5_results_Pl$best_pars$fL,
                                  dH = model5_results_Pl$best_pars$dH, eH = model5_results_Pl$best_pars$eH, fH = model5_results_Pl$best_pars$fH,
                                  Xo = model5_results_Pl$best_pars$Xo, Xb = model5_results_Pl$best_pars$Xb,
                                  SNR = x, temp = mean(model5_results_Pl$source_data$DD5), aSMR = mean(model5_results_Pl$source_data$aSMR))
lines(x,y)


#MaxGrowth predictions from model. how well do they fit SBSmc2?
SI_conversion$predMaxGrowth <-rep(NA, 4)
SI_conversion$predMaxGrowth[1]<-SORTIE.climate.sens(MaxGrowth=model5_results_Pl$best_pars$MaxGrowth , aL = model5_results_Pl$best_pars$aL, bL = model5_results_Pl$best_pars$bL, cL = model5_results_Pl$best_pars$cL,
                    aH = model5_results_Pl$best_pars$aH, bH = model5_results_Pl$best_pars$bH, cH = model5_results_Pl$best_pars$cH,
                    dL = model5_results_Pl$best_pars$dL, eL = model5_results_Pl$best_pars$eL, fL = model5_results_Pl$best_pars$fL,
                    dH = model5_results_Pl$best_pars$dH, eH = model5_results_Pl$best_pars$eH, fH = model5_results_Pl$best_pars$fH,
                    Xo = model5_results_Pl$best_pars$Xo, Xb = model5_results_Pl$best_pars$Xb,
                    SNR = 2, temp = 875.73, aSMR = 3)
SI_conversion$predMaxGrowth[2]<-SORTIE.climate.sens(MaxGrowth=model5_results_Pl$best_pars$MaxGrowth , aL = model5_results_Pl$best_pars$aL, bL = model5_results_Pl$best_pars$bL, cL = model5_results_Pl$best_pars$cL,
                                                    aH = model5_results_Pl$best_pars$aH, bH = model5_results_Pl$best_pars$bH, cH = model5_results_Pl$best_pars$cH,
                                                    dL = model5_results_Pl$best_pars$dL, eL = model5_results_Pl$best_pars$eL, fL = model5_results_Pl$best_pars$fL,
                                                    dH = model5_results_Pl$best_pars$dH, eH = model5_results_Pl$best_pars$eH, fH = model5_results_Pl$best_pars$fH,
                                                    Xo = model5_results_Pl$best_pars$Xo, Xb = model5_results_Pl$best_pars$Xb,
                                                    SNR = 3, temp = 875.73, aSMR = 4)
SI_conversion$predMaxGrowth[3]<-SORTIE.climate.sens(MaxGrowth=model5_results_Pl$best_pars$MaxGrowth , aL = model5_results_Pl$best_pars$aL, bL = model5_results_Pl$best_pars$bL, cL = model5_results_Pl$best_pars$cL,
                                                    aH = model5_results_Pl$best_pars$aH, bH = model5_results_Pl$best_pars$bH, cH = model5_results_Pl$best_pars$cH,
                                                    dL = model5_results_Pl$best_pars$dL, eL = model5_results_Pl$best_pars$eL, fL = model5_results_Pl$best_pars$fL,
                                                    dH = model5_results_Pl$best_pars$dH, eH = model5_results_Pl$best_pars$eH, fH = model5_results_Pl$best_pars$fH,
                                                    Xo = model5_results_Pl$best_pars$Xo, Xb = model5_results_Pl$best_pars$Xb,
                                                    SNR = 4, temp = 875.73, aSMR = 5)
SI_conversion$predMaxGrowth[4]<-SORTIE.climate.sens(MaxGrowth=model5_results_Pl$best_pars$MaxGrowth , aL = model5_results_Pl$best_pars$aL, bL = model5_results_Pl$best_pars$bL, cL = model5_results_Pl$best_pars$cL,
                                                    aH = model5_results_Pl$best_pars$aH, bH = model5_results_Pl$best_pars$bH, cH = model5_results_Pl$best_pars$cH,
                                                    dL = model5_results_Pl$best_pars$dL, eL = model5_results_Pl$best_pars$eL, fL = model5_results_Pl$best_pars$fL,
                                                    dH = model5_results_Pl$best_pars$dH, eH = model5_results_Pl$best_pars$eH, fH = model5_results_Pl$best_pars$fH,
                                                    Xo = model5_results_Pl$best_pars$Xo, Xb = model5_results_Pl$best_pars$Xb,
                                                    SNR = 5, temp = 875.73, aSMR = 6)
SI_conversion


###################################################################################################################
###################################################################################################################
####################################        Subalpine fir           ######################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

#######################
############################# SORTIE aSMR double logistic; double logistic for temp and including SNR
######################  fitting without noESSFwk1   which gives it a big jog
######################
SI_climate_Bl.noESSFwk1<-subset(SI_climate_Bl, SI_climate_Bl$BGC !=  "ESSFwk1")

SORTIE.aSMRSNR.temp.model.log2 <- function(maxSI, aL, bL, cL, aH, bH, cH, aSMR, dL, eL, fL, dH, eH, fH,  temp, Xo, Xb, SNR)
{  maxSI * (aL+ (1-aL)/(1+ ((bL/aSMR)^cL))) * (aH+ (1-aH)/(1+ ((aSMR/bH)^cH)) ) * (dL+ (1-dL)/(1+ ((eL/temp)^fL))) * (dH+ (1-dH)/(1+ ((temp/eH)^fH)) ) * exp(-0.5 * ((SNR - Xo)/Xb)^2) }

var <- list(aSMR = "aSMR", temp = "DD5", SNR = "SNR" )
par <- list( maxSI = 42, aL = 0.001, bL = 3, cL = .5, aH =.05, bH = 8, cH = 70, dL = 0.1, eL = 500, fL = 1, dH =4, eH = 500, fH = 0.01,  Xo = 3, Xb = 1, sd = 1)
par_lo <- list( maxSI = 0,  aL = 0,     bL = 0,  cL = 0, aH =0, bH = 0,  cH = 0, dL = 0,  eL = 0,   fL = 0,  dH =0,  eH = 0,    fH = 0, Xo = 0, Xb = 0,    sd = 0)
par_hi <- list( maxSI = 150, aL = .01, bL = 6, cL = 10, aH =.2, bH = 15, cH = 100, dL = 1, eL = 2000, fL = 20, dH =1, eH = 2000, fH = 500, Xo = 50, Xb = 40,   sd = 10)

var$x <- "SI" #this is the site index
var$mean <- "predicted"
var$log <- TRUE

model_results_Bl <- anneal (model = SORTIE.aSMRSNR.temp.model.log2, par = par, var = var, source_data = SI_climate_Bl.noESSFwk1, par_lo = par_lo, par_hi = par_hi,
                            pdf=dnorm, dep_var="SI",   initial_temp = 5, temp_red = 0.96, max_iter=100000, hessian = FALSE)

write_results(model_results_Bl,"SORTIE parameterization/Bl1_SiteIndex.txt", data=F)

#look at residuals
plot(model_results_Bl$source_data$predicted, model_results_Bl$source_data$SI)
SBSmc2_points<-subset(model_results_Bl$source_data, model_results_Bl$source_data$BGC == "SBSmc2")
points(SBSmc2_points$predicted, SBSmc2_points$SI, pch = 2, col = "red")
abline(0,1)
xyplot(predicted~ SI|BGC, data = model_results_Bl$source_data)
xyplot(predicted~ SI, group = BGC, data = model_results_Bl$source_data)

resid<- model_results_Bl$source_data$predicted-model_results_Bl$source_data$SI
plot(model_results_Bl$source_data$SI, resid)
##  plot a 0 line  (i.e. intercept = 0, slope = 0)
abline(0,0)

plot(model_results_Bl$source_data$aSMR, model_results_Bl$source_data$SI)
x<-seq(1, 8, .1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Bl$best_pars$maxSI , aL = model_results_Bl$best_pars$aL, bL = model_results_Bl$best_pars$bL, cL = model_results_Bl$best_pars$cL,
                                  aH = model_results_Bl$best_pars$aH, bH = model_results_Bl$best_pars$bH, cH = model_results_Bl$best_pars$cH,
                                  dL = model_results_Bl$best_pars$dL, eL = model_results_Bl$best_pars$eL, fL = model_results_Bl$best_pars$fL,
                                  dH = model_results_Bl$best_pars$dH, eH = model_results_Bl$best_pars$eH, fH = model_results_Bl$best_pars$fH,
                                  Xo = model_results_Bl$best_pars$Xo, Xb = model_results_Bl$best_pars$Xb,
                                  SNR = mean(model_results_Bl$source_data$SNR), aSMR = x, temp = mean(model_results_Bl$source_data$DD5))
lines(x,y)

plot(model_results_Bl$source_data$DD5, model_results_Bl$source_data$SI, xlim = c(300, 1500))
x<-seq(300, 1500, 1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Bl$best_pars$maxSI , aL = model_results_Bl$best_pars$aL, bL = model_results_Bl$best_pars$bL, cL = model_results_Bl$best_pars$cL,
                                  aH = model_results_Bl$best_pars$aH, bH = model_results_Bl$best_pars$bH, cH = model_results_Bl$best_pars$cH,
                                  dL = model_results_Bl$best_pars$dL, eL = model_results_Bl$best_pars$eL, fL = model_results_Bl$best_pars$fL,
                                  dH = model_results_Bl$best_pars$dH, eH = model_results_Bl$best_pars$eH, fH = model_results_Bl$best_pars$fH,
                                  Xo = model_results_Bl$best_pars$Xo, Xb = model_results_Bl$best_pars$Xb,
                                  SNR = mean(model_results_Bl$source_data$SNR), temp = x, aSMR = mean(model_results_Bl$source_data$aSMR))
lines(x,y)
SBSmc2_points<-subset(SI_climate_Bl, SI_climate_Bl$BGC == "SBSmc2")
points(SBSmc2_points$DD5, SBSmc2_points$SI, pch = 2, col = "red")
SI_climate_Bl[which(SI_climate_Bl$DD5 >1150),]

plot(model_results_Bl$source_data$SNR, model_results_Bl$source_data$SI)
x<-seq(1, 5, .1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Bl$best_pars$maxSI , aL = model_results_Bl$best_pars$aL, bL = model_results_Bl$best_pars$bL, cL = model_results_Bl$best_pars$cL,
                                  aH = model_results_Bl$best_pars$aH, bH = model_results_Bl$best_pars$bH, cH = model_results_Bl$best_pars$cH,
                                  dL = model_results_Bl$best_pars$dL, eL = model_results_Bl$best_pars$eL, fL = model_results_Bl$best_pars$fL,
                                  dH = model_results_Bl$best_pars$dH, eH = model_results_Bl$best_pars$eH, fH = model_results_Bl$best_pars$fH,
                                  Xo = model_results_Bl$best_pars$Xo, Xb = model_results_Bl$best_pars$Xb,
                                  SNR = x, temp = mean(model_results_Bl$source_data$DD5), aSMR = mean(model_results_Bl$source_data$aSMR))
lines(x,y)
SBSmc2_points<-subset(SI_climate_Bl, SI_climate_Bl$BGC == "SBSmc2")
points(SBSmc2_points$SNR, SBSmc2_points$SI, pch = 2, col = "red")

###################
####################
####################
#take above function and predict new MaxGrowth parameters
###################
####################
####################
subset(SMRCross, SMRCross$BGC == "SBSmc2")

#Create conversion from Site Index to MaxGrowth
SI_conversion<-matrix(nrow = 4,ncol = 2)

#predict site index for SBSmc2 '09' sites...
SI_conversion[4,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Bl$best_pars$maxSI , aL = model_results_Bl$best_pars$aL, bL = model_results_Bl$best_pars$bL, cL = model_results_Bl$best_pars$cL,
                                                   aH = model_results_Bl$best_pars$aH, bH = model_results_Bl$best_pars$bH, cH = model_results_Bl$best_pars$cH,
                                                   dL = model_results_Bl$best_pars$dL, eL = model_results_Bl$best_pars$eL, fL = model_results_Bl$best_pars$fL,
                                                   dH = model_results_Bl$best_pars$dH, eH = model_results_Bl$best_pars$eH, fH = model_results_Bl$best_pars$fH,
                                                   Xo = model_results_Bl$best_pars$Xo, Xb = model_results_Bl$best_pars$Xb,
                                                   SNR = 5, temp = 875.73, aSMR = 6)
#predict site index for SBSmc2 '06' sites
SI_conversion[3,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Bl$best_pars$maxSI , aL = model_results_Bl$best_pars$aL, bL = model_results_Bl$best_pars$bL, cL = model_results_Bl$best_pars$cL,
                                                   aH = model_results_Bl$best_pars$aH, bH = model_results_Bl$best_pars$bH, cH = model_results_Bl$best_pars$cH,
                                                   dL = model_results_Bl$best_pars$dL, eL = model_results_Bl$best_pars$eL, fL = model_results_Bl$best_pars$fL,
                                                   dH = model_results_Bl$best_pars$dH, eH = model_results_Bl$best_pars$eH, fH = model_results_Bl$best_pars$fH,
                                                   Xo = model_results_Bl$best_pars$Xo, Xb = model_results_Bl$best_pars$Xb,
                                                   SNR = 4, temp = 875.73, aSMR = 5)
#predict site index for SBSmc2 '01' sites
SI_conversion[2,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Bl$best_pars$maxSI , aL = model_results_Bl$best_pars$aL, bL = model_results_Bl$best_pars$bL, cL = model_results_Bl$best_pars$cL,
                                                   aH = model_results_Bl$best_pars$aH, bH = model_results_Bl$best_pars$bH, cH = model_results_Bl$best_pars$cH,
                                                   dL = model_results_Bl$best_pars$dL, eL = model_results_Bl$best_pars$eL, fL = model_results_Bl$best_pars$fL,
                                                   dH = model_results_Bl$best_pars$dH, eH = model_results_Bl$best_pars$eH, fH = model_results_Bl$best_pars$fH,
                                                   Xo = model_results_Bl$best_pars$Xo, Xb = model_results_Bl$best_pars$Xb,
                                                   SNR = 3, temp = 875.73, aSMR = 4)

#predict site index for SBSmc2 '02' sites
SI_conversion[1,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Bl$best_pars$maxSI , aL = model_results_Bl$best_pars$aL, bL = model_results_Bl$best_pars$bL, cL = model_results_Bl$best_pars$cL,
                                                   aH = model_results_Bl$best_pars$aH, bH = model_results_Bl$best_pars$bH, cH = model_results_Bl$best_pars$cH,
                                                   dL = model_results_Bl$best_pars$dL, eL = model_results_Bl$best_pars$eL, fL = model_results_Bl$best_pars$fL,
                                                   dH = model_results_Bl$best_pars$dH, eH = model_results_Bl$best_pars$eH, fH = model_results_Bl$best_pars$fH,
                                                   Xo = model_results_Bl$best_pars$Xo, Xb = model_results_Bl$best_pars$Xb,
                                                   SNR = 2, temp = 875.73, aSMR = 3)

#Adult MaxGrowth for subalpine fir is 4.3 (mm/yr) * (SNR/6)^0.65
SI_conversion[1,2]<- 4.3 * (2/5)^0.65
SI_conversion[2,2]<- 4.3 * (3/5)^0.65
SI_conversion[3,2]<- 4.3 * (4/5)^0.65
SI_conversion[4,2]<- 4.3 * (5/5)^0.65
colnames(SI_conversion)<- c("SI", "MaxGrowth")
SI_conversion<-as.data.frame(SI_conversion)
plot(SI_conversion$SI, SI_conversion$MaxGrowth)
slope<- (SI_conversion$MaxGrowth[4]-SI_conversion$MaxGrowth[1])/(SI_conversion$SI[4]-SI_conversion$SI[1])
# y = mx + b
SI_conversion$MaxGrowth[4] = SI_conversion$SI[4] * slope + b
SI_conversion$MaxGrowth[4]-SI_conversion$SI[4] * slope
slope
SI_conv_fnct<-function(x) {0.8205226 * x + -9.88971}
x<- seq(13, 21, .1)
y<- SI_conv_fnct(x)
lines(x, y)

SI_climate_Bl.noESSFwk1$MaxGrowth_conv<- SI_conv_fnct(SI_climate_Bl.noESSFwk1$SI)

##################                                ####################
################## relative productivity fitting ###################
#################                                 ###################

SORTIE.climate.sens <- function(MaxGrowth, aL, bL, cL, aH, bH, cH, aSMR, dL, eL, fL, dH, eH, fH,  temp, Xo, Xb, SNR)
{  MaxGrowth * (aL+ (1-aL)/(1+ ((bL/aSMR)^cL))) * (aH+ (1-aH)/(1+ ((aSMR/bH)^cH)) ) * (dL+ (1-dL)/(1+ ((eL/temp)^fL))) * (dH+ (1-dH)/(1+ ((temp/eH)^fH)) ) * exp(-0.5 * ((SNR - Xo)/Xb)^2) }

var <- list(aSMR = "aSMR", temp = "DD5", SNR = "SNR" )
par <- list( MaxGrowth = 1, aL = 0.001, bL = 3, cL = .5, aH =.05, bH = 8, cH = 70, dL = 0.1, eL = 500, fL = 1, dH =4, eH = 500, fH = 0.01,  Xo = 3, Xb = 1, sd = 1)
par_lo <- list( MaxGrowth = 0,  aL = 0,     bL = 0,  cL = 0, aH =0, bH = 0,  cH = 0, dL = 0,  eL = 0,   fL = 0,  dH =0,  eH = 0,    fH = 0, Xo = 0, Xb = 0,    sd = 0)
par_hi <- list( MaxGrowth = 60, aL = .05, bL = 6, cL = 10, aH =.1, bH = 15, cH = 50, dL = 1, eL = 2000, fL = 30, dH =1, eH = 2000, fH = 500, Xo = 40, Xb = 35,   sd = 10)

var$x <- "MaxGrowth_conv" #this is the site index
var$mean <- "predicted"
var$log <- TRUE

model5_results_Bl <- anneal (model = SORTIE.climate.sens, par = par, var = var, source_data = SI_climate_Bl.noESSFwk1, par_lo = par_lo, par_hi = par_hi,
                             pdf=dnorm, dep_var="MaxGrowth_conv",   initial_temp = 5, temp_red = 0.96, max_iter=100000, hessian = FALSE)

write_results(model5_results_Bl,"SORTIE parameterization/Bl3_MaxGrowth_adult.txt", data=F)

#look at residuals
plot(model5_results_Bl$source_data$predicted, model5_results_Bl$source_data$MaxGrowth_conv)
SBSmc2_points<-subset(model5_results_Bl$source_data, model5_results_Bl$source_data$BGC == "SBSmc2")
points(SBSmc2_points$predicted, SBSmc2_points$MaxGrowth_conv, pch = 2, col = "red")
abline(0,1)
xyplot(predicted~ MaxGrowth_conv|BGC, data = model5_results_Bl$source_data)
xyplot(predicted~ MaxGrowth_conv, group = BGC, data = model5_results_Bl$source_data)

resid<- model5_results_Bl$source_data$predicted-model5_results_Bl$source_data$MaxGrowth_conv
plot(model5_results_Bl$source_data$MaxGrowth_conv, resid)
##  plot a 0 line  (i.e. intercept = 0, slope = 0)
abline(0,0)

plot(model5_results_Bl$source_data$aSMR, model5_results_Bl$source_data$MaxGrowth_conv)
x<-seq(1, 8, .1)
y<-SORTIE.climate.sens(MaxGrowth=model5_results_Bl$best_pars$MaxGrowth , aL = model5_results_Bl$best_pars$aL, bL = model5_results_Bl$best_pars$bL, cL = model5_results_Bl$best_pars$cL,
                       aH = model5_results_Bl$best_pars$aH, bH = model5_results_Bl$best_pars$bH, cH = model5_results_Bl$best_pars$cH,
                       dL = model5_results_Bl$best_pars$dL, eL = model5_results_Bl$best_pars$eL, fL = model5_results_Bl$best_pars$fL,
                       dH = model5_results_Bl$best_pars$dH, eH = model5_results_Bl$best_pars$eH, fH = model5_results_Bl$best_pars$fH,
                       Xo = model5_results_Bl$best_pars$Xo, Xb = model5_results_Bl$best_pars$Xb,
                       SNR = mean(model5_results_Bl$source_data$SNR), aSMR = x, temp = mean(model5_results_Bl$source_data$DD5))
lines(x,y)

plot(model5_results_Bl$source_data$DD5, model5_results_Bl$source_data$MaxGrowth_conv)
x<-seq(300, 1500, 1)
y<-SORTIE.climate.sens(MaxGrowth=model5_results_Bl$best_pars$MaxGrowth , aL = model5_results_Bl$best_pars$aL, bL = model5_results_Bl$best_pars$bL, cL = model5_results_Bl$best_pars$cL,
                       aH = model5_results_Bl$best_pars$aH, bH = model5_results_Bl$best_pars$bH, cH = model5_results_Bl$best_pars$cH,
                       dL = model5_results_Bl$best_pars$dL, eL = model5_results_Bl$best_pars$eL, fL = model5_results_Bl$best_pars$fL,
                       dH = model5_results_Bl$best_pars$dH, eH = model5_results_Bl$best_pars$eH, fH = model5_results_Bl$best_pars$fH,
                       Xo = model5_results_Bl$best_pars$Xo, Xb = model5_results_Bl$best_pars$Xb,
                       SNR = mean(model5_results_Bl$source_data$SNR), temp = x, aSMR = mean(model5_results_Bl$source_data$aSMR))
lines(x,y)
#can a new maxgrowth be put in site index function? Answer is No
SI_conv_fnct(model5_results_Bl$best_pars$MaxGrowth)
y2<-SORTIE.aSMRSNR.temp.model.log2(maxSI=9.07  , aL = model_results_Bl$best_pars$aL, bL = model_results_Bl$best_pars$bL, cL = model_results_Bl$best_pars$cL,
                                   aH = model_results_Bl$best_pars$aH, bH = model_results_Bl$best_pars$bH, cH = model_results_Bl$best_pars$cH,
                                   dL = model_results_Bl$best_pars$dL, eL = model_results_Bl$best_pars$eL, fL = model_results_Bl$best_pars$fL,
                                   dH = model_results_Bl$best_pars$dH, eH = model_results_Bl$best_pars$eH, fH = model_results_Bl$best_pars$fH,
                                   Xo = model_results_Bl$best_pars$Xo, Xb = model_results_Bl$best_pars$Xb,
                                   SNR = mean(model5_results_Bl$source_data$SNR), temp = x, aSMR = mean(model5_results_Bl$source_data$aSMR))
lines(x,y2, lwd = 2, col = "red")


plot(model5_results_Bl$source_data$SNR, model5_results_Bl$source_data$MaxGrowth_conv)
x<-seq(1, 5, .1)
y<-SORTIE.climate.sens(MaxGrowth=model5_results_Bl$best_pars$MaxGrowth , aL = model5_results_Bl$best_pars$aL, bL = model5_results_Bl$best_pars$bL, cL = model5_results_Bl$best_pars$cL,
                       aH = model5_results_Bl$best_pars$aH, bH = model5_results_Bl$best_pars$bH, cH = model5_results_Bl$best_pars$cH,
                       dL = model5_results_Bl$best_pars$dL, eL = model5_results_Bl$best_pars$eL, fL = model5_results_Bl$best_pars$fL,
                       dH = model5_results_Bl$best_pars$dH, eH = model5_results_Bl$best_pars$eH, fH = model5_results_Bl$best_pars$fH,
                       Xo = model5_results_Bl$best_pars$Xo, Xb = model5_results_Bl$best_pars$Xb,
                       SNR = x, temp = mean(model5_results_Bl$source_data$DD5), aSMR = mean(model5_results_Bl$source_data$aSMR))
lines(x,y)


#MaxGrowth predictions from model. how well do they fit SBSmc2?
SI_conversion$predMaxGrowth <-rep(NA, 4)
SI_conversion$predMaxGrowth[1]<-SORTIE.climate.sens(MaxGrowth=model5_results_Bl$best_pars$MaxGrowth , aL = model5_results_Bl$best_pars$aL, bL = model5_results_Bl$best_pars$bL, cL = model5_results_Bl$best_pars$cL,
                                                    aH = model5_results_Bl$best_pars$aH, bH = model5_results_Bl$best_pars$bH, cH = model5_results_Bl$best_pars$cH,
                                                    dL = model5_results_Bl$best_pars$dL, eL = model5_results_Bl$best_pars$eL, fL = model5_results_Bl$best_pars$fL,
                                                    dH = model5_results_Bl$best_pars$dH, eH = model5_results_Bl$best_pars$eH, fH = model5_results_Bl$best_pars$fH,
                                                    Xo = model5_results_Bl$best_pars$Xo, Xb = model5_results_Bl$best_pars$Xb,
                                                    SNR = 2, temp = 875.73, aSMR = 3)
SI_conversion$predMaxGrowth[2]<-SORTIE.climate.sens(MaxGrowth=model5_results_Bl$best_pars$MaxGrowth , aL = model5_results_Bl$best_pars$aL, bL = model5_results_Bl$best_pars$bL, cL = model5_results_Bl$best_pars$cL,
                                                    aH = model5_results_Bl$best_pars$aH, bH = model5_results_Bl$best_pars$bH, cH = model5_results_Bl$best_pars$cH,
                                                    dL = model5_results_Bl$best_pars$dL, eL = model5_results_Bl$best_pars$eL, fL = model5_results_Bl$best_pars$fL,
                                                    dH = model5_results_Bl$best_pars$dH, eH = model5_results_Bl$best_pars$eH, fH = model5_results_Bl$best_pars$fH,
                                                    Xo = model5_results_Bl$best_pars$Xo, Xb = model5_results_Bl$best_pars$Xb,
                                                    SNR = 3, temp = 875.73, aSMR = 4)
SI_conversion$predMaxGrowth[3]<-SORTIE.climate.sens(MaxGrowth=model5_results_Bl$best_pars$MaxGrowth , aL = model5_results_Bl$best_pars$aL, bL = model5_results_Bl$best_pars$bL, cL = model5_results_Bl$best_pars$cL,
                                                    aH = model5_results_Bl$best_pars$aH, bH = model5_results_Bl$best_pars$bH, cH = model5_results_Bl$best_pars$cH,
                                                    dL = model5_results_Bl$best_pars$dL, eL = model5_results_Bl$best_pars$eL, fL = model5_results_Bl$best_pars$fL,
                                                    dH = model5_results_Bl$best_pars$dH, eH = model5_results_Bl$best_pars$eH, fH = model5_results_Bl$best_pars$fH,
                                                    Xo = model5_results_Bl$best_pars$Xo, Xb = model5_results_Bl$best_pars$Xb,
                                                    SNR = 4, temp = 875.73, aSMR = 5)
SI_conversion$predMaxGrowth[4]<-SORTIE.climate.sens(MaxGrowth=model5_results_Bl$best_pars$MaxGrowth , aL = model5_results_Bl$best_pars$aL, bL = model5_results_Bl$best_pars$bL, cL = model5_results_Bl$best_pars$cL,
                                                    aH = model5_results_Bl$best_pars$aH, bH = model5_results_Bl$best_pars$bH, cH = model5_results_Bl$best_pars$cH,
                                                    dL = model5_results_Bl$best_pars$dL, eL = model5_results_Bl$best_pars$eL, fL = model5_results_Bl$best_pars$fL,
                                                    dH = model5_results_Bl$best_pars$dH, eH = model5_results_Bl$best_pars$eH, fH = model5_results_Bl$best_pars$fH,
                                                    Xo = model5_results_Bl$best_pars$Xo, Xb = model5_results_Bl$best_pars$Xb,
                                                    SNR = 5, temp = 875.73, aSMR = 6)
SI_conversion


###################################################################################################################
###################################################################################################################
####################################         Spruce           ######################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

#######################
############################# SORTIE aSMR double logistic; double logistic for temp and including SNR
######################
######################
SI_climate_Sx.noICHwk4<-subset(SI_climate_Sx, SI_climate_Sx$BGC !=  "ICHwk4")

SORTIE.aSMRSNR.temp.model.log2 <- function(maxSI, aL, bL, cL, aH, bH, cH, aSMR, dL, eL, fL, dH, eH, fH,  temp, Xo, Xb, SNR)
{  maxSI * (aL+ (1-aL)/(1+ ((bL/aSMR)^cL))) * (aH+ (1-aH)/(1+ ((aSMR/bH)^cH)) ) * (dL+ (1-dL)/(1+ ((eL/temp)^fL))) * (dH+ (1-dH)/(1+ ((temp/eH)^fH)) ) * exp(-0.5 * ((SNR - Xo)/Xb)^2) }

var <- list(aSMR = "aSMR", temp = "DD5", SNR = "SNR" )
par <- list( maxSI = 42, aL = 0.001, bL = 3, cL = .5, aH =.05, bH = 8, cH = 70, dL = 0.1, eL = 500, fL = 1, dH =4, eH = 500, fH = 0.01,  Xo = 3, Xb = 1, sd = 1)
par_lo <- list( maxSI = 0,  aL = 0,     bL = 0,  cL = 0, aH =0, bH = 0,  cH = 0, dL = 0,  eL = 0,   fL = 0,  dH =0,  eH = 0,    fH = 0, Xo = 0, Xb = 0,    sd = 0)
par_hi <- list( maxSI = 150, aL = .01, bL = 6, cL = 10, aH =.1, bH = 15, cH = 100, dL = 1, eL = 3000, fL = 20, dH =1, eH = 2000, fH = 500, Xo = 50, Xb = 35,   sd = 10)

var$x <- "SI" #this is the site index
var$mean <- "predicted"
var$log <- TRUE

model_results_Sx <- anneal (model = SORTIE.aSMRSNR.temp.model.log2, par = par, var = var, source_data = SI_climate_Sx.noICHwk4, par_lo = par_lo, par_hi = par_hi,
                            pdf=dnorm, dep_var="SI",   initial_temp = 5, temp_red = 0.96, max_iter=100000, hessian = FALSE)

write_results(model_results_Sx,"SORTIE parameterization/Sx1_SiteIndex.txt", data=F)

#look at residuals
plot(model_results_Sx$source_data$predicted, model_results_Sx$source_data$SI)
SBSmc2_points<-subset(model_results_Sx$source_data, model_results_Sx$source_data$BGC == "SBSmc2")
points(SBSmc2_points$predicted, SBSmc2_points$SI, pch = 2, col = "red")
abline(0,1)
xyplot(predicted~ SI|BGC, data = model_results_Sx$source_data)
xyplot(predicted~ SI, group = BGC, data = model_results_Sx$source_data)

resid<- model_results_Sx$source_data$predicted-model_results_Sx$source_data$SI
plot(model_results_Sx$source_data$SI, resid)
##  plot a 0 line  (i.e. intercept = 0, slope = 0)
abline(0,0)

plot(model_results_Sx$source_data$aSMR, model_results_Sx$source_data$SI)
x<-seq(1, 8, .1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Sx$best_pars$maxSI , aL = model_results_Sx$best_pars$aL, bL = model_results_Sx$best_pars$bL, cL = model_results_Sx$best_pars$cL,
                                  aH = model_results_Sx$best_pars$aH, bH = model_results_Sx$best_pars$bH, cH = model_results_Sx$best_pars$cH,
                                  dL = model_results_Sx$best_pars$dL, eL = model_results_Sx$best_pars$eL, fL = model_results_Sx$best_pars$fL,
                                  dH = model_results_Sx$best_pars$dH, eH = model_results_Sx$best_pars$eH, fH = model_results_Sx$best_pars$fH,
                                  Xo = model_results_Sx$best_pars$Xo, Xb = model_results_Sx$best_pars$Xb,
                                  SNR = mean(model_results_Sx$source_data$SNR), aSMR = x, temp = mean(model_results_Sx$source_data$DD5))
lines(x,y)

plot(model_results_Sx$source_data$DD5, model_results_Sx$source_data$SI, xlim = c(300, 1500))
x<-seq(300, 1500, 1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Sx$best_pars$maxSI , aL = model_results_Sx$best_pars$aL, bL = model_results_Sx$best_pars$bL, cL = model_results_Sx$best_pars$cL,
                                  aH = model_results_Sx$best_pars$aH, bH = model_results_Sx$best_pars$bH, cH = model_results_Sx$best_pars$cH,
                                  dL = model_results_Sx$best_pars$dL, eL = model_results_Sx$best_pars$eL, fL = model_results_Sx$best_pars$fL,
                                  dH = model_results_Sx$best_pars$dH, eH = model_results_Sx$best_pars$eH, fH = model_results_Sx$best_pars$fH,
                                  Xo = model_results_Sx$best_pars$Xo, Xb = model_results_Sx$best_pars$Xb,
                                  SNR = mean(model_results_Sx$source_data$SNR), temp = x, aSMR = mean(model_results_Sx$source_data$aSMR))
lines(x,y)
SBSmc2_points<-subset(SI_climate_Sx, SI_climate_Sx$BGC == "SBSmc2")
points(SBSmc2_points$DD5, SBSmc2_points$SI, pch = 2, col = "red")
SI_climate_Sx[which(SI_climate_Sx$DD5 >1050),]

plot(model_results_Sx$source_data$SNR, model_results_Sx$source_data$SI)
x<-seq(1, 5, .1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Sx$best_pars$maxSI , aL = model_results_Sx$best_pars$aL, bL = model_results_Sx$best_pars$bL, cL = model_results_Sx$best_pars$cL,
                                  aH = model_results_Sx$best_pars$aH, bH = model_results_Sx$best_pars$bH, cH = model_results_Sx$best_pars$cH,
                                  dL = model_results_Sx$best_pars$dL, eL = model_results_Sx$best_pars$eL, fL = model_results_Sx$best_pars$fL,
                                  dH = model_results_Sx$best_pars$dH, eH = model_results_Sx$best_pars$eH, fH = model_results_Sx$best_pars$fH,
                                  Xo = model_results_Sx$best_pars$Xo, Xb = model_results_Sx$best_pars$Xb,
                                  SNR = x, temp = mean(model_results_Sx$source_data$DD5), aSMR = mean(model_results_Sx$source_data$aSMR))
lines(x,y)
SBSmc2_points<-subset(SI_climate_Sx, SI_climate_Sx$BGC == "SBSmc2")
points(SBSmc2_points$SNR, SBSmc2_points$SI, pch = 2, col = "red")


 ###################
####################
####################
#take above function and predict new MaxGrowth parameters
###################
####################
####################
subset(SMRCross, SMRCross$BGC == "SBSmc2")

#Create conversion from Site Index to MaxGrowth
SI_conversion<-matrix(nrow = 4,ncol = 2)

#predict site index for SBSmc2 '09' sites...
SI_conversion[4,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Sx$best_pars$maxSI , aL = model_results_Sx$best_pars$aL, bL = model_results_Sx$best_pars$bL, cL = model_results_Sx$best_pars$cL,
                                                   aH = model_results_Sx$best_pars$aH, bH = model_results_Sx$best_pars$bH, cH = model_results_Sx$best_pars$cH,
                                                   dL = model_results_Sx$best_pars$dL, eL = model_results_Sx$best_pars$eL, fL = model_results_Sx$best_pars$fL,
                                                   dH = model_results_Sx$best_pars$dH, eH = model_results_Sx$best_pars$eH, fH = model_results_Sx$best_pars$fH,
                                                   Xo = model_results_Sx$best_pars$Xo, Xb = model_results_Sx$best_pars$Xb,
                                                   SNR = 5, temp = 875.73, aSMR = 6)
#predict site index for SBSmc2 '06' sites
SI_conversion[3,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Sx$best_pars$maxSI , aL = model_results_Sx$best_pars$aL, bL = model_results_Sx$best_pars$bL, cL = model_results_Sx$best_pars$cL,
                                                   aH = model_results_Sx$best_pars$aH, bH = model_results_Sx$best_pars$bH, cH = model_results_Sx$best_pars$cH,
                                                   dL = model_results_Sx$best_pars$dL, eL = model_results_Sx$best_pars$eL, fL = model_results_Sx$best_pars$fL,
                                                   dH = model_results_Sx$best_pars$dH, eH = model_results_Sx$best_pars$eH, fH = model_results_Sx$best_pars$fH,
                                                   Xo = model_results_Sx$best_pars$Xo, Xb = model_results_Sx$best_pars$Xb,
                                                   SNR = 4, temp = 875.73, aSMR = 5)
#predict site index for SBSmc2 '01' sites
SI_conversion[2,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Sx$best_pars$maxSI , aL = model_results_Sx$best_pars$aL, bL = model_results_Sx$best_pars$bL, cL = model_results_Sx$best_pars$cL,
                                                   aH = model_results_Sx$best_pars$aH, bH = model_results_Sx$best_pars$bH, cH = model_results_Sx$best_pars$cH,
                                                   dL = model_results_Sx$best_pars$dL, eL = model_results_Sx$best_pars$eL, fL = model_results_Sx$best_pars$fL,
                                                   dH = model_results_Sx$best_pars$dH, eH = model_results_Sx$best_pars$eH, fH = model_results_Sx$best_pars$fH,
                                                   Xo = model_results_Sx$best_pars$Xo, Xb = model_results_Sx$best_pars$Xb,
                                                   SNR = 3, temp = 875.73, aSMR = 4)

#predict site index for SBSmc2 '02' sites
SI_conversion[1,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Sx$best_pars$maxSI , aL = model_results_Sx$best_pars$aL, bL = model_results_Sx$best_pars$bL, cL = model_results_Sx$best_pars$cL,
                                                   aH = model_results_Sx$best_pars$aH, bH = model_results_Sx$best_pars$bH, cH = model_results_Sx$best_pars$cH,
                                                   dL = model_results_Sx$best_pars$dL, eL = model_results_Sx$best_pars$eL, fL = model_results_Sx$best_pars$fL,
                                                   dH = model_results_Sx$best_pars$dH, eH = model_results_Sx$best_pars$eH, fH = model_results_Sx$best_pars$fH,
                                                   Xo = model_results_Sx$best_pars$Xo, Xb = model_results_Sx$best_pars$Xb,
                                                   SNR = 2, temp = 875.73, aSMR = 3)

#Adult MaxGrowth for spruce is 7.7 (mm/yr) * (SNR/6)^1.5
SI_conversion[1,2]<- 7.7 * (2/5)^1.5
SI_conversion[2,2]<- 7.7 * (3/5)^1.5
SI_conversion[3,2]<- 7.7 * (4/5)^1.5
SI_conversion[4,2]<- 7.7 * (5/5)^1.5
colnames(SI_conversion)<- c("SI", "MaxGrowth")
SI_conversion<-as.data.frame(SI_conversion)
plot(SI_conversion$SI, SI_conversion$MaxGrowth)
slope<- (SI_conversion$MaxGrowth[4]-SI_conversion$MaxGrowth[1])/(SI_conversion$SI[4]-SI_conversion$SI[1])
# y = mx + b
SI_conversion$MaxGrowth[4] = SI_conversion$SI[4] * slope + b
SI_conversion$MaxGrowth[4]-SI_conversion$SI[4] * slope
slope
#first conversion - SI_conv_fnct<-function(x) {1.556549 * x + -22.2432}
SI_conv_fnct<-function(x) {1.499451 * x + -21.19854} #this is the 4th conversion
x<- seq(13, 21, .1)
y<- SI_conv_fnct(x)
lines(x, y)

SI_climate_Sx.noICHwk4$MaxGrowth_conv<- SI_conv_fnct(SI_climate_Sx.noICHwk4$SI)

#add a couple more points for conversion fitting non-linear
SI_conversion[5,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Sx$best_pars$maxSI , aL = model_results_Sx$best_pars$aL, bL = model_results_Sx$best_pars$bL, cL = model_results_Sx$best_pars$cL,
                                                   aH = model_results_Sx$best_pars$aH, bH = model_results_Sx$best_pars$bH, cH = model_results_Sx$best_pars$cH,
                                                   dL = model_results_Sx$best_pars$dL, eL = model_results_Sx$best_pars$eL, fL = model_results_Sx$best_pars$fL,
                                                   dH = model_results_Sx$best_pars$dH, eH = model_results_Sx$best_pars$eH, fH = model_results_Sx$best_pars$fH,
                                                   Xo = model_results_Sx$best_pars$Xo, Xb = model_results_Sx$best_pars$Xb,
                                                   SNR = 1.5, temp = 875.73, aSMR = 2.5)
SI_conversion[5,2]<- 7.7 * (1.5/5)^1.5

SI_conversion[6,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Sx$best_pars$maxSI , aL = model_results_Sx$best_pars$aL, bL = model_results_Sx$best_pars$bL, cL = model_results_Sx$best_pars$cL,
                                                   aH = model_results_Sx$best_pars$aH, bH = model_results_Sx$best_pars$bH, cH = model_results_Sx$best_pars$cH,
                                                   dL = model_results_Sx$best_pars$dL, eL = model_results_Sx$best_pars$eL, fL = model_results_Sx$best_pars$fL,
                                                   dH = model_results_Sx$best_pars$dH, eH = model_results_Sx$best_pars$eH, fH = model_results_Sx$best_pars$fH,
                                                   Xo = model_results_Sx$best_pars$Xo, Xb = model_results_Sx$best_pars$Xb,
                                                   SNR = 4.5, temp = 875.73, aSMR = 5.5)
SI_conversion[6,2]<- 7.7 * (4.5/5)^1.5


#try conversion without 09 site because with 09 site predictions aren't great
slope<- (SI_conversion$MaxGrowth[3]-SI_conversion$MaxGrowth[1])/(SI_conversion$SI[3]-SI_conversion$SI[1])
SI_conversion$MaxGrowth[3]-SI_conversion$SI[3] * slope
slope
SI_conv_fnct<-function(x) {1.315794 * x + -18.50149}
x<- seq(13, 21, .1)
y2<- SI_conv_fnct(x)
lines(x, y2, lwd = 2)

SI_climate_Sx.noICHwk4$MaxGrowth_conv<- SI_conv_fnct(SI_climate_Sx.noICHwk4$SI)


#try non-linear conversion
#add a couple more points for conversion fitting non-linear
SI_conversion[5,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Sx$best_pars$maxSI , aL = model_results_Sx$best_pars$aL, bL = model_results_Sx$best_pars$bL, cL = model_results_Sx$best_pars$cL,
                                                   aH = model_results_Sx$best_pars$aH, bH = model_results_Sx$best_pars$bH, cH = model_results_Sx$best_pars$cH,
                                                   dL = model_results_Sx$best_pars$dL, eL = model_results_Sx$best_pars$eL, fL = model_results_Sx$best_pars$fL,
                                                   dH = model_results_Sx$best_pars$dH, eH = model_results_Sx$best_pars$eH, fH = model_results_Sx$best_pars$fH,
                                                   Xo = model_results_Sx$best_pars$Xo, Xb = model_results_Sx$best_pars$Xb,
                                                   SNR = 1.5, temp = 875.73, aSMR = 2.5)
SI_conversion[5,2]<- 7.7 * (1.5/5)^1.5

SI_conversion[6,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_Sx$best_pars$maxSI , aL = model_results_Sx$best_pars$aL, bL = model_results_Sx$best_pars$bL, cL = model_results_Sx$best_pars$cL,
                                                   aH = model_results_Sx$best_pars$aH, bH = model_results_Sx$best_pars$bH, cH = model_results_Sx$best_pars$cH,
                                                   dL = model_results_Sx$best_pars$dL, eL = model_results_Sx$best_pars$eL, fL = model_results_Sx$best_pars$fL,
                                                   dH = model_results_Sx$best_pars$dH, eH = model_results_Sx$best_pars$eH, fH = model_results_Sx$best_pars$fH,
                                                   Xo = model_results_Sx$best_pars$Xo, Xb = model_results_Sx$best_pars$Xb,
                                                   SNR = 4.5, temp = 875.73, aSMR = 5.5)
SI_conversion[6,2]<- 7.7 * (4.5/5)^1.5

plot(SI_conversion$SI, SI_conversion$MaxGrowth, ylim = c(0,8), xlim = c(0,20))


##############
############## fitting non-linear model
#############

conv.fit <- function(b, SI, c) {  b + exp(SI*c) }

var <- list( SI = "SI" )
par <- list( b = -20, c = 1)
par_lo <- list(  b = -40, c = -2)
par_hi <- list(  b = 40, c = 2)

var$x <- "MaxGrowth" #this is the site index
var$mean <- "predicted"
var$log <- TRUE

conv.fit_results_Sx <- anneal (model = conv.fit, par = par, var = var, source_data = SI_conversion, par_lo = par_lo, par_hi = par_hi,
                               pdf=dnorm, dep_var="MaxGrowth",   initial_temp = 5, temp_red = 0.90, max_iter=10000, hessian = FALSE)


plot(SI_conversion$SI, SI_conversion$MaxGrowth)
x<- seq(13, 21, .1)
y<- conv.fit(b = conv.fit_results_Sx$best_pars$b, c = conv.fit_results_Sx$best_pars$c, SI = x)
lines(x, y)

SI_conv_fnct<-function(x) { exp(x*0.136237) + -6.278732} #3rd conversion
SI_climate_Sx.noICHwk4$MaxGrowth_conv<- SI_conv_fnct(SI_climate_Sx.noICHwk4$SI)
SI_conv_fnct<-function(x) { exp(x*0.1344173) + -5.951995} #5th conversion


##################                                ####################
################## relative productivity fitting ###################
#################                                 ###################

SORTIE.climate.sens <- function(MaxGrowth, aL, bL, cL, aH, bH, cH, aSMR, dL, eL, fL, dH, eH, fH,  temp, Xo, Xb, SNR)
{  MaxGrowth * (aL+ (1-aL)/(1+ ((bL/aSMR)^cL))) * (aH+ (1-aH)/(1+ ((aSMR/bH)^cH)) ) * (dL+ (1-dL)/(1+ ((eL/temp)^fL))) * (dH+ (1-dH)/(1+ ((temp/eH)^fH)) ) * exp(-0.5 * ((SNR - Xo)/Xb)^2) }

var <- list(aSMR = "aSMR", temp = "DD5", SNR = "SNR" )
par <- list( MaxGrowth = 40, aL = 0.001, bL = 3, cL = .5, aH =.05, bH = 8, cH = 70, dL = 0.1, eL = 500, fL = 1, dH =4, eH = 500, fH = 0.01,  Xo = 3, Xb = 1, sd = 1)
par_lo <- list( MaxGrowth = 0,  aL = 0,     bL = 0,  cL = 0, aH =0, bH = 0,  cH = 0, dL = 0,  eL = 0,   fL = 0,  dH =0,  eH = 0,    fH = 0, Xo = 0, Xb = 0,    sd = 0)
par_hi <- list( MaxGrowth = 200, aL = .01, bL = 6, cL = 10, aH =.1, bH = 15, cH = 50, dL = 1, eL = 2000, fL = 20, dH =1, eH = 2000, fH = 500, Xo = 40, Xb = 35,   sd = 10)

var$x <- "MaxGrowth_conv" #this is the site index
var$mean <- "predicted"
var$log <- TRUE

model5_results_Sx <- anneal (model = SORTIE.climate.sens, par = par, var = var, source_data = SI_climate_Sx.noICHwk4, par_lo = par_lo, par_hi = par_hi,
                             pdf=dnorm, dep_var="MaxGrowth_conv",   initial_temp = 5, temp_red = 0.96, max_iter=100000, hessian = FALSE)

write_results(model5_results_Sx,"SORTIE parameterization/Sx3_MaxGrowth_adult_5th_conv.txt", data=F)

#look at residuals
plot(model5_results_Sx$source_data$predicted, model5_results_Sx$source_data$MaxGrowth_conv)
SBSmc2_points<-subset(model5_results_Sx$source_data, model5_results_Sx$source_data$BGC == "SBSmc2")
points(SBSmc2_points$predicted, SBSmc2_points$MaxGrowth_conv, pch = 2, col = "red")
abline(0,1)
xyplot(predicted~ MaxGrowth_conv|BGC, data = model5_results_Sx$source_data)
xyplot(predicted~ MaxGrowth_conv, group = BGC, data = model5_results_Sx$source_data)

resid<- model5_results_Sx$source_data$predicted-model5_results_Sx$source_data$MaxGrowth_conv
plot(model5_results_Sx$source_data$MaxGrowth_conv, resid)
##  plot a 0 line  (i.e. intercept = 0, slope = 0)
abline(0,0)

plot(model5_results_Sx$source_data$aSMR, model5_results_Sx$source_data$MaxGrowth_conv)
x<-seq(1, 8, .1)
y<-SORTIE.climate.sens(MaxGrowth=model5_results_Sx$best_pars$MaxGrowth , aL = model5_results_Sx$best_pars$aL, bL = model5_results_Sx$best_pars$bL, cL = model5_results_Sx$best_pars$cL,
                       aH = model5_results_Sx$best_pars$aH, bH = model5_results_Sx$best_pars$bH, cH = model5_results_Sx$best_pars$cH,
                       dL = model5_results_Sx$best_pars$dL, eL = model5_results_Sx$best_pars$eL, fL = model5_results_Sx$best_pars$fL,
                       dH = model5_results_Sx$best_pars$dH, eH = model5_results_Sx$best_pars$eH, fH = model5_results_Sx$best_pars$fH,
                       Xo = model5_results_Sx$best_pars$Xo, Xb = model5_results_Sx$best_pars$Xb,
                       SNR = mean(model5_results_Sx$source_data$SNR), aSMR = x, temp = mean(model5_results_Sx$source_data$DD5))
lines(x,y)

plot(model5_results_Sx$source_data$DD5, model5_results_Sx$source_data$MaxGrowth_conv, xlim = c(300, 1500))
x<-seq(300, 1500, 1)
y<-SORTIE.climate.sens(MaxGrowth=model5_results_Sx$best_pars$MaxGrowth , aL = model5_results_Sx$best_pars$aL, bL = model5_results_Sx$best_pars$bL, cL = model5_results_Sx$best_pars$cL,
                       aH = model5_results_Sx$best_pars$aH, bH = model5_results_Sx$best_pars$bH, cH = model5_results_Sx$best_pars$cH,
                       dL = model5_results_Sx$best_pars$dL, eL = model5_results_Sx$best_pars$eL, fL = model5_results_Sx$best_pars$fL,
                       dH = model5_results_Sx$best_pars$dH, eH = model5_results_Sx$best_pars$eH, fH = model5_results_Sx$best_pars$fH,
                       Xo = model5_results_Sx$best_pars$Xo, Xb = model5_results_Sx$best_pars$Xb,
                       SNR = mean(model5_results_Sx$source_data$SNR), temp = x, aSMR = mean(model5_results_Sx$source_data$aSMR))
lines(x,y)

plot(model5_results_Sx$source_data$SNR, model5_results_Sx$source_data$MaxGrowth_conv)
x<-seq(1, 5, .1)
y<-SORTIE.climate.sens(MaxGrowth=model5_results_Sx$best_pars$MaxGrowth , aL = model5_results_Sx$best_pars$aL, bL = model5_results_Sx$best_pars$bL, cL = model5_results_Sx$best_pars$cL,
                       aH = model5_results_Sx$best_pars$aH, bH = model5_results_Sx$best_pars$bH, cH = model5_results_Sx$best_pars$cH,
                       dL = model5_results_Sx$best_pars$dL, eL = model5_results_Sx$best_pars$eL, fL = model5_results_Sx$best_pars$fL,
                       dH = model5_results_Sx$best_pars$dH, eH = model5_results_Sx$best_pars$eH, fH = model5_results_Sx$best_pars$fH,
                       Xo = model5_results_Sx$best_pars$Xo, Xb = model5_results_Sx$best_pars$Xb,
                       SNR = x, temp = mean(model5_results_Sx$source_data$DD5), aSMR = mean(model5_results_Sx$source_data$aSMR))
lines(x,y)


#MaxGrowth predictions from model. how well do they fit SBSmc2?
SI_conversion$predMaxGrowth <-rep(NA, 6)
SI_conversion$predMaxGrowth[1]<-SORTIE.climate.sens(MaxGrowth=model5_results_Sx$best_pars$MaxGrowth , aL = model5_results_Sx$best_pars$aL, bL = model5_results_Sx$best_pars$bL, cL = model5_results_Sx$best_pars$cL,
                                                    aH = model5_results_Sx$best_pars$aH, bH = model5_results_Sx$best_pars$bH, cH = model5_results_Sx$best_pars$cH,
                                                    dL = model5_results_Sx$best_pars$dL, eL = model5_results_Sx$best_pars$eL, fL = model5_results_Sx$best_pars$fL,
                                                    dH = model5_results_Sx$best_pars$dH, eH = model5_results_Sx$best_pars$eH, fH = model5_results_Sx$best_pars$fH,
                                                    Xo = model5_results_Sx$best_pars$Xo, Xb = model5_results_Sx$best_pars$Xb,
                                                    SNR = 2, temp = 875.73, aSMR = 3)
SI_conversion$predMaxGrowth[2]<-SORTIE.climate.sens(MaxGrowth=model5_results_Sx$best_pars$MaxGrowth , aL = model5_results_Sx$best_pars$aL, bL = model5_results_Sx$best_pars$bL, cL = model5_results_Sx$best_pars$cL,
                                                    aH = model5_results_Sx$best_pars$aH, bH = model5_results_Sx$best_pars$bH, cH = model5_results_Sx$best_pars$cH,
                                                    dL = model5_results_Sx$best_pars$dL, eL = model5_results_Sx$best_pars$eL, fL = model5_results_Sx$best_pars$fL,
                                                    dH = model5_results_Sx$best_pars$dH, eH = model5_results_Sx$best_pars$eH, fH = model5_results_Sx$best_pars$fH,
                                                    Xo = model5_results_Sx$best_pars$Xo, Xb = model5_results_Sx$best_pars$Xb,
                                                    SNR = 3, temp = 875.73, aSMR = 4)
SI_conversion$predMaxGrowth[3]<-SORTIE.climate.sens(MaxGrowth=model5_results_Sx$best_pars$MaxGrowth , aL = model5_results_Sx$best_pars$aL, bL = model5_results_Sx$best_pars$bL, cL = model5_results_Sx$best_pars$cL,
                                                    aH = model5_results_Sx$best_pars$aH, bH = model5_results_Sx$best_pars$bH, cH = model5_results_Sx$best_pars$cH,
                                                    dL = model5_results_Sx$best_pars$dL, eL = model5_results_Sx$best_pars$eL, fL = model5_results_Sx$best_pars$fL,
                                                    dH = model5_results_Sx$best_pars$dH, eH = model5_results_Sx$best_pars$eH, fH = model5_results_Sx$best_pars$fH,
                                                    Xo = model5_results_Sx$best_pars$Xo, Xb = model5_results_Sx$best_pars$Xb,
                                                    SNR = 4, temp = 875.73, aSMR = 5)
SI_conversion$predMaxGrowth[4]<-SORTIE.climate.sens(MaxGrowth=model5_results_Sx$best_pars$MaxGrowth , aL = model5_results_Sx$best_pars$aL, bL = model5_results_Sx$best_pars$bL, cL = model5_results_Sx$best_pars$cL,
                                                    aH = model5_results_Sx$best_pars$aH, bH = model5_results_Sx$best_pars$bH, cH = model5_results_Sx$best_pars$cH,
                                                    dL = model5_results_Sx$best_pars$dL, eL = model5_results_Sx$best_pars$eL, fL = model5_results_Sx$best_pars$fL,
                                                    dH = model5_results_Sx$best_pars$dH, eH = model5_results_Sx$best_pars$eH, fH = model5_results_Sx$best_pars$fH,
                                                    Xo = model5_results_Sx$best_pars$Xo, Xb = model5_results_Sx$best_pars$Xb,
                                                    SNR = 5, temp = 875.73, aSMR = 6)
SI_conversion

plot(SI_conversion$MaxGrowth, SI_conversion$predMaxGrowth)
abline(0,1)


###################################################################################################################
###################################################################################################################
####################################         Aspen           ######################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

#######################
############################# SORTIE aSMR double logistic; double logistic for temp and including SNR
######################
######################

SORTIE.aSMRSNR.temp.model.log2 <- function(maxSI, aL, bL, cL, aH, bH, cH, aSMR, dL, eL, fL, dH, eH, fH,  temp, Xo, Xb, SNR)
{  maxSI * (aL+ (1-aL)/(1+ ((bL/aSMR)^cL))) * (aH+ (1-aH)/(1+ ((aSMR/bH)^cH)) ) * (dL+ (1-dL)/(1+ ((eL/temp)^fL))) * (dH+ (1-dH)/(1+ ((temp/eH)^fH)) ) * exp(-0.5 * ((SNR - Xo)/Xb)^2) }

var <- list(aSMR = "aSMR", temp = "DD5", SNR = "SNR" )
par <- list( maxSI = 42, aL = 0.001, bL = 3, cL = .5, aH =.05, bH = 8, cH = 70, dL = 0.1, eL = 500, fL = 1, dH =4, eH = 500, fH = 0.01,  Xo = 3, Xb = 1, sd = 1)
par_lo <- list( maxSI = 0,  aL = 0,     bL = 0,  cL = 0, aH =0, bH = 0,  cH = 0, dL = 0,  eL = 0,   fL = 0,  dH =0,  eH = 0,    fH = 0, Xo = 0, Xb = 0,    sd = 0)
par_hi <- list( maxSI = 150, aL = .01, bL = 6, cL = 10, aH =.1, bH = 15, cH = 100, dL = 1, eL = 2000, fL = 20, dH =1, eH = 2000, fH = 500, Xo = 50, Xb = 45,   sd = 10)

var$x <- "SI" #this is the site index
var$mean <- "predicted"
var$log <- TRUE

model_results_At <- anneal (model = SORTIE.aSMRSNR.temp.model.log2, par = par, var = var, source_data = SI_climate_At, par_lo = par_lo, par_hi = par_hi,
                            pdf=dnorm, dep_var="SI",   initial_temp = 5, temp_red = 0.96, max_iter=100000, hessian = FALSE)

write_results(model_results_At,"SORTIE parameterization/At1_SiteIndex.txt", data=F)

#look at residuals
plot(model_results_At$source_data$predicted, model_results_At$source_data$SI)
SBSmc2_points<-subset(model_results_At$source_data, model_results_At$source_data$BGC == "SBSmc2")
points(SBSmc2_points$predicted, SBSmc2_points$SI, pch = 2, col = "red")
abline(0,1)
xyplot(predicted~ SI|BGC, data = model_results_At$source_data)
xyplot(predicted~ SI, group = BGC, data = model_results_At$source_data)

resid<- model_results_At$source_data$predicted-model_results_At$source_data$SI
plot(model_results_At$source_data$SI, resid)
##  plot a 0 line  (i.e. intercept = 0, slope = 0)
abline(0,0)

plot(model_results_At$source_data$aSMR, model_results_At$source_data$SI)
x<-seq(1, 8, .1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_At$best_pars$maxSI , aL = model_results_At$best_pars$aL, bL = model_results_At$best_pars$bL, cL = model_results_At$best_pars$cL,
                                  aH = model_results_At$best_pars$aH, bH = model_results_At$best_pars$bH, cH = model_results_At$best_pars$cH,
                                  dL = model_results_At$best_pars$dL, eL = model_results_At$best_pars$eL, fL = model_results_At$best_pars$fL,
                                  dH = model_results_At$best_pars$dH, eH = model_results_At$best_pars$eH, fH = model_results_At$best_pars$fH,
                                  Xo = model_results_At$best_pars$Xo, Xb = model_results_At$best_pars$Xb,
                                  SNR = mean(model_results_At$source_data$SNR), aSMR = x, temp = mean(model_results_At$source_data$DD5))
lines(x,y)

plot(model_results_At$source_data$DD5, model_results_At$source_data$SI, xlim = c(300, 1500))
x<-seq(300, 1500, 1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_At$best_pars$maxSI , aL = model_results_At$best_pars$aL, bL = model_results_At$best_pars$bL, cL = model_results_At$best_pars$cL,
                                  aH = model_results_At$best_pars$aH, bH = model_results_At$best_pars$bH, cH = model_results_At$best_pars$cH,
                                  dL = model_results_At$best_pars$dL, eL = model_results_At$best_pars$eL, fL = model_results_At$best_pars$fL,
                                  dH = model_results_At$best_pars$dH, eH = model_results_At$best_pars$eH, fH = model_results_At$best_pars$fH,
                                  Xo = model_results_At$best_pars$Xo, Xb = model_results_At$best_pars$Xb,
                                  SNR = mean(model_results_At$source_data$SNR), temp = x, aSMR = mean(model_results_At$source_data$aSMR))
lines(x,y)

plot(model_results_At$source_data$SNR, model_results_At$source_data$SI)
x<-seq(1, 5, .1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_At$best_pars$maxSI , aL = model_results_At$best_pars$aL, bL = model_results_At$best_pars$bL, cL = model_results_At$best_pars$cL,
                                  aH = model_results_At$best_pars$aH, bH = model_results_At$best_pars$bH, cH = model_results_At$best_pars$cH,
                                  dL = model_results_At$best_pars$dL, eL = model_results_At$best_pars$eL, fL = model_results_At$best_pars$fL,
                                  dH = model_results_At$best_pars$dH, eH = model_results_At$best_pars$eH, fH = model_results_At$best_pars$fH,
                                  Xo = model_results_At$best_pars$Xo, Xb = model_results_At$best_pars$Xb,
                                  SNR = x, temp = mean(model_results_At$source_data$DD5), aSMR = mean(model_results_At$source_data$aSMR))
lines(x,y)
SBSmc2_points<-subset(SI_climate_At, SI_climate_At$BGC == "SBSmc2")
points(SBSmc2_points$SNR, SBSmc2_points$SI, pch = 2, col = "red")


#create fake data from Chen et al. 2002 to fill in a wider DD5
#Model 2 from chen et al. 2002
#SI = BWBS × (29.5 – 0.019(ELE) – 2.12(P) + 3.12(F_SMR) + 6.21(M_SMR) – 7.36(RG)
#             – 3.56(N) – 2.67(W)) + SB × (17.85 – 2.8(FL) + 4.6(M_SMR)) + IDF
#× (340.95 – 2.694(LON) + 2.89(SD_SMR) – 6.95(S)) + MS × (24.3 – 4.32(MD_SMR)
#+ 3.94(R) + 2.74(N) – 5.29(S)) + ICH × (12.6 + 0.011(ELE) + 4.61(M_SMR) + 5.37(FL))
# ELE,elevation; LAT, latitude; LON, longitude. Soil moisture regime: MD_SMR, moderately dry; SD_SMR, slightly dry; F_SMR, fresh; and M_SMR, moist.
#Soil nutrient regime: VP, very poor and P, poor. Slope–aspect: RG, ridge; E, east; S, south; W, west; N, north; and FL, flat.

BWBS_chen<-function(ELE, P, F_SMR, M_SMR, RG, N, W) {(29.5 - 0.019*(ELE) - 2.12*(P) + 3.12*(F_SMR) + 6.21*(M_SMR) - 7.36*(RG)    - 3.56*(N) - 2.67*(W))}
SB_chen<-function(FL, M_SMR) {(17.85 - 2.8*(FL) + 4.6*(M_SMR))}
IDF_chen<- function(LON, SD_SMR, S) {(340.95 - 2.694*(LON) + 2.89*(SD_SMR) - 6.95*(S))    }
MS_chen<-function(MD_SMR, R, N, S) {24.3 - 4.32*(MD_SMR)+ 3.94*(R) + 2.74*(N) - 5.29*(S)}
ICH_chen<- function(ELE, M_SMR, FL) {(12.6 + 0.011*(ELE) + 4.61*(M_SMR) + 5.37*(FL))}

#example spot is unit B4 at Date Creek
ICH_test_aSMR6<-ICH_chen(ELE = 498, M_SMR = 1, FL =0)
ICH_test_aSMR6
ICH_test_aSMR5<-ICH_chen(ELE = 498, M_SMR = 0, FL =0)
ICH_test_aSMR5

head(SI_climate_At)
new_SI_climate_At<-SI_climate_At[1,]
new_SI_climate_At[1,] <-c("ICH_chen", NA, 4, 6, 22.688, 498, 1208.72, 209.8, NA)
new_SI_climate_At[2,] <-c("ICH_chen", NA, 3, 5, 18.078, 498, 1208.72, 209.8, NA)

#example SB spot is smithers community forest, logging road at nordic centre
SB_test_aSMR6<-SB_chen( M_SMR = 1, FL =0)
SB_test_aSMR6
SB_test_aSMR5<-SB_chen( M_SMR = 0, FL =0)
SB_test_aSMR5
new_SI_climate_At[3,] <-c("SB_chen", NA, 4, 6, 22.45, 910, 875.73, 201.7, NA)
new_SI_climate_At[4,] <-c("SB_chen", NA, 3, 5, 17.85, 910, 875.73, 201.7, NA)

#example BWBS spot is kinaskan lake.. pretend west aspect
BWBS_test_aSMR6<-BWBS_chen( ELE =818, P =0, F_SMR =0, M_SMR =1, RG =0, N =0, W =1)
BWBS_test_aSMR6
BWBS_test_aSMR5<-BWBS_chen( ELE =818, P =0, F_SMR =1, M_SMR =0, RG =0, N =0, W =1)
BWBS_test_aSMR5
BWBS_test_aSMR4<-BWBS_chen( ELE =818, P =0, F_SMR =0, M_SMR =0, RG =0, N =0, W =1)
BWBS_test_aSMR4
BWBS_test_aSMR3<-BWBS_chen( ELE =818, P =1, F_SMR =0, M_SMR =0, RG =0, N =0, W =1)
BWBS_test_aSMR3
subset(SMRCross, SMRCross$BGC == "BWBSdk")
new_SI_climate_At[5,] <-c("BWBS_chen", NA, 4, 6, 17.498, 818, 696.94, 158.47, NA)
new_SI_climate_At[6,] <-c("BWBS_chen", NA, 3, 5, 14.408, 818, 696.94, 158.47, NA)
new_SI_climate_At[7,] <-c("BWBS_chen", NA, 3, 4, 11.288, 818, 696.94, 158.47, NA)
new_SI_climate_At[8,] <-c("BWBS_chen", NA, 2, 4, 9.168, 818, 696.94, 158.47, NA)

#example IDF = IDFdk3 by williams lake airport
subset(SMRCross, SMRCross$BGC == "IDFdk3")
IDF_chen<- function(LON, SD_SMR, S) {(340.95 - 2.694*(LON) + 2.89*(SD_SMR) - 6.95*(S))    }
IDF_test_aSMR2<-IDF_chen( LON = 122, SD_SMR = 1, S = 0)
IDF_test_aSMR2
new_SI_climate_At[9,] <-c("IDF_chen", NA, 2, 2, 15.172, 939, 1113.71, 295.44, NA)

#example MSdw - 49.86278	-115.89474	1367	MSdw near Kimberley
subset(SMRCross, SMRCross$BGC == "MSdw")
MS_chen<-function(MD_SMR, R, N, S) {24.3 - 4.32*(MD_SMR)+ 3.94*(R) + 2.74*(N) - 5.29*(S)}
MS_test_aSMR3<-MS_chen( MD_SMR = 1, R = 0, N = 0, S = 0)
MS_test_aSMR3
MS_test_aSMR4<-MS_chen( MD_SMR = 0, R = 0, N = 0, S = 0)
MS_test_aSMR4
new_SI_climate_At[10,] <-c("MS_chen", NA, 2, 3, 19.98, 1367, 1130.81, 215.42, NA)
new_SI_climate_At[11,] <-c("MS_chen", NA, 3, 4, 24.3, 1367, 1130.81, 215.42, NA)

new_SI_climate_At$SNR<-as.numeric(new_SI_climate_At$SNR)
new_SI_climate_At$X<-as.numeric(new_SI_climate_At$X)
new_SI_climate_At$aSMR<-as.numeric(new_SI_climate_At$aSMR)
new_SI_climate_At$SI<-as.numeric(new_SI_climate_At$SI)
new_SI_climate_At$Elevation<-as.numeric(new_SI_climate_At$Elevation)
new_SI_climate_At$DD5<-as.numeric(new_SI_climate_At$DD5)
new_SI_climate_At$CMD<-as.numeric(new_SI_climate_At$CMD)
new_SI_climate_At$MaxGrowth_conv<-as.numeric(new_SI_climate_At$MaxGrowth_conv)
SI_climate_Atplus<-rbind(SI_climate_At, new_SI_climate_At)

model_results_At_new <- anneal (model = SORTIE.aSMRSNR.temp.model.log2, par = par, var = var, source_data = SI_climate_Atplus, par_lo = par_lo, par_hi = par_hi,
                            pdf=dnorm, dep_var="SI",   initial_temp = 5, temp_red = 0.96, max_iter=100000, hessian = FALSE)

write_results(model_results_At_new,"SORTIE parameterization/At1new_SiteIndex.txt", data=F)

#look at residuals
plot(model_results_At_new$source_data$predicted, model_results_At_new$source_data$SI)
SBSmc2_points<-subset(model_results_At_new$source_data, model_results_At_new$source_data$BGC == "SBSmc2")
points(SBSmc2_points$predicted, SBSmc2_points$SI, pch = 2, col = "red")
abline(0,1)

plot(model_results_At_new$source_data$aSMR, model_results_At_new$source_data$SI)
x<-seq(1, 8, .1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_At_new$best_pars$maxSI , aL = model_results_At_new$best_pars$aL, bL = model_results_At_new$best_pars$bL, cL = model_results_At_new$best_pars$cL,
                                  aH = model_results_At_new$best_pars$aH, bH = model_results_At_new$best_pars$bH, cH = model_results_At_new$best_pars$cH,
                                  dL = model_results_At_new$best_pars$dL, eL = model_results_At_new$best_pars$eL, fL = model_results_At_new$best_pars$fL,
                                  dH = model_results_At_new$best_pars$dH, eH = model_results_At_new$best_pars$eH, fH = model_results_At_new$best_pars$fH,
                                  Xo = model_results_At_new$best_pars$Xo, Xb = model_results_At_new$best_pars$Xb,
                                  SNR = mean(model_results_At_new$source_data$SNR), aSMR = x, temp = mean(model_results_At_new$source_data$DD5))
lines(x,y)

plot(model_results_At_new$source_data$DD5, model_results_At_new$source_data$SI, xlim = c(300, 1500))
x<-seq(300, 1500, 1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_At_new$best_pars$maxSI , aL = model_results_At_new$best_pars$aL, bL = model_results_At_new$best_pars$bL, cL = model_results_At_new$best_pars$cL,
                                  aH = model_results_At_new$best_pars$aH, bH = model_results_At_new$best_pars$bH, cH = model_results_At_new$best_pars$cH,
                                  dL = model_results_At_new$best_pars$dL, eL = model_results_At_new$best_pars$eL, fL = model_results_At_new$best_pars$fL,
                                  dH = model_results_At_new$best_pars$dH, eH = model_results_At_new$best_pars$eH, fH = model_results_At_new$best_pars$fH,
                                  Xo = model_results_At_new$best_pars$Xo, Xb = model_results_At_new$best_pars$Xb,
                                  SNR = mean(model_results_At_new$source_data$SNR), temp = x, aSMR = mean(model_results_At_new$source_data$aSMR))
lines(x,y)

plot(model_results_At_new$source_data$SNR, model_results_At_new$source_data$SI)
x<-seq(1, 5, .1)
y<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_At_new$best_pars$maxSI , aL = model_results_At_new$best_pars$aL, bL = model_results_At_new$best_pars$bL, cL = model_results_At_new$best_pars$cL,
                                  aH = model_results_At_new$best_pars$aH, bH = model_results_At_new$best_pars$bH, cH = model_results_At_new$best_pars$cH,
                                  dL = model_results_At_new$best_pars$dL, eL = model_results_At_new$best_pars$eL, fL = model_results_At_new$best_pars$fL,
                                  dH = model_results_At_new$best_pars$dH, eH = model_results_At_new$best_pars$eH, fH = model_results_At_new$best_pars$fH,
                                  Xo = model_results_At_new$best_pars$Xo, Xb = model_results_At_new$best_pars$Xb,
                                  SNR = x, temp = mean(model_results_At_new$source_data$DD5), aSMR = mean(model_results_At_new$source_data$aSMR))
lines(x,y)

###################
####################
####################
#take above function and predict new MaxGrowth parameters
###################
####################
####################
subset(SMRCross, SMRCross$BGC == "SBSmc2")

#Create conversion from Site Index to MaxGrowth
SI_conversion<-matrix(nrow = 4,ncol = 2)

#predict site index for SBSmc2 '09' sites...
SI_conversion[4,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_At_new$best_pars$maxSI , aL = model_results_At_new$best_pars$aL, bL = model_results_At_new$best_pars$bL, cL = model_results_At_new$best_pars$cL,
                                                   aH = model_results_At_new$best_pars$aH, bH = model_results_At_new$best_pars$bH, cH = model_results_At_new$best_pars$cH,
                                                   dL = model_results_At_new$best_pars$dL, eL = model_results_At_new$best_pars$eL, fL = model_results_At_new$best_pars$fL,
                                                   dH = model_results_At_new$best_pars$dH, eH = model_results_At_new$best_pars$eH, fH = model_results_At_new$best_pars$fH,
                                                   Xo = model_results_At_new$best_pars$Xo, Xb = model_results_At_new$best_pars$Xb,
                                                   SNR = 5, temp = 875.73, aSMR = 6)
#predict site index for SBSmc2 '06' sites
SI_conversion[3,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_At_new$best_pars$maxSI , aL = model_results_At_new$best_pars$aL, bL = model_results_At_new$best_pars$bL, cL = model_results_At_new$best_pars$cL,
                                                   aH = model_results_At_new$best_pars$aH, bH = model_results_At_new$best_pars$bH, cH = model_results_At_new$best_pars$cH,
                                                   dL = model_results_At_new$best_pars$dL, eL = model_results_At_new$best_pars$eL, fL = model_results_At_new$best_pars$fL,
                                                   dH = model_results_At_new$best_pars$dH, eH = model_results_At_new$best_pars$eH, fH = model_results_At_new$best_pars$fH,
                                                   Xo = model_results_At_new$best_pars$Xo, Xb = model_results_At_new$best_pars$Xb,
                                                   SNR = 4, temp = 875.73, aSMR = 5)
#predict site index for SBSmc2 '01' sites
SI_conversion[2,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_At_new$best_pars$maxSI , aL = model_results_At_new$best_pars$aL, bL = model_results_At_new$best_pars$bL, cL = model_results_At_new$best_pars$cL,
                                                   aH = model_results_At_new$best_pars$aH, bH = model_results_At_new$best_pars$bH, cH = model_results_At_new$best_pars$cH,
                                                   dL = model_results_At_new$best_pars$dL, eL = model_results_At_new$best_pars$eL, fL = model_results_At_new$best_pars$fL,
                                                   dH = model_results_At_new$best_pars$dH, eH = model_results_At_new$best_pars$eH, fH = model_results_At_new$best_pars$fH,
                                                   Xo = model_results_At_new$best_pars$Xo, Xb = model_results_At_new$best_pars$Xb,
                                                   SNR = 3, temp = 875.73, aSMR = 4)

#predict site index for SBSmc2 '02' sites
SI_conversion[1,1]<-SORTIE.aSMRSNR.temp.model.log2(maxSI=model_results_At_new$best_pars$maxSI , aL = model_results_At_new$best_pars$aL, bL = model_results_At_new$best_pars$bL, cL = model_results_At_new$best_pars$cL,
                                                   aH = model_results_At_new$best_pars$aH, bH = model_results_At_new$best_pars$bH, cH = model_results_At_new$best_pars$cH,
                                                   dL = model_results_At_new$best_pars$dL, eL = model_results_At_new$best_pars$eL, fL = model_results_At_new$best_pars$fL,
                                                   dH = model_results_At_new$best_pars$dH, eH = model_results_At_new$best_pars$eH, fH = model_results_At_new$best_pars$fH,
                                                   Xo = model_results_At_new$best_pars$Xo, Xb = model_results_At_new$best_pars$Xb,
                                                   SNR = 2, temp = 875.73, aSMR = 3)

#Adult MaxGrowth for aspen is 4.8 (mm/yr) * (SNR/6)^1.1
SI_conversion[1,2]<- 4.8 * (2/5)^1.1
SI_conversion[2,2]<- 4.8 * (3/5)^1.1
SI_conversion[3,2]<- 4.8 * (4/5)^1.1
SI_conversion[4,2]<- 4.8 * (5/5)^1.1
colnames(SI_conversion)<- c("SI", "MaxGrowth")
SI_conversion<-as.data.frame(SI_conversion)
plot(SI_conversion$SI, SI_conversion$MaxGrowth)
slope<- (SI_conversion$MaxGrowth[4]-SI_conversion$MaxGrowth[1])/(SI_conversion$SI[4]-SI_conversion$SI[1])
# y = mx + b
SI_conversion$MaxGrowth[4] = SI_conversion$SI[4] * slope + b
SI_conversion$MaxGrowth[4]-SI_conversion$SI[4] * slope
slope
SI_conv_fnct<-function(x) {0.6164963 * x + -7.319682}
x<- seq(13, 21, .1)
y<- SI_conv_fnct(x)
lines(x, y)

SI_climate_At$MaxGrowth_conv<- SI_conv_fnct(SI_climate_At$SI)
SI_climate_Atplus$MaxGrowth_conv<- SI_conv_fnct(SI_climate_Atplus$SI)

##################                                ####################
################## relative productivity fitting ###################
#################                                 ###################

SORTIE.climate.sens <- function(MaxGrowth, aL, bL, cL, aH, bH, cH, aSMR, dL, eL, fL, dH, eH, fH,  temp, Xo, Xb, SNR)
{  MaxGrowth * (aL+ (1-aL)/(1+ ((bL/aSMR)^cL))) * (aH+ (1-aH)/(1+ ((aSMR/bH)^cH)) ) * (dL+ (1-dL)/(1+ ((eL/temp)^fL))) * (dH+ (1-dH)/(1+ ((temp/eH)^fH)) ) * exp(-0.5 * ((SNR - Xo)/Xb)^2) }

var <- list(aSMR = "aSMR", temp = "DD5", SNR = "SNR" )
par <- list( MaxGrowth = 1, aL = 0.001, bL = 3, cL = .5, aH =.05, bH = 8, cH = 70, dL = 0.1, eL = 500, fL = 1, dH =4, eH = 500, fH = 0.01,  Xo = 3, Xb = 1, sd = 1)
par_lo <- list( MaxGrowth = 0,  aL = 0,     bL = 0,  cL = 0, aH =0, bH = 0,  cH = 0, dL = 0,  eL = 0,   fL = 0,  dH =0,  eH = 0,    fH = 0, Xo = 0, Xb = 0,    sd = 0)
par_hi <- list( MaxGrowth = 60, aL = .01, bL = 6, cL = 10, aH =.1, bH = 15, cH = 50, dL = 1, eL = 2000, fL = 20, dH =1, eH = 2000, fH = 500, Xo = 50, Xb = 35,   sd = 10)

var$x <- "MaxGrowth_conv" #this is the site index
var$mean <- "predicted"
var$log <- TRUE

model5_results_At <- anneal (model = SORTIE.climate.sens, par = par, var = var, source_data = SI_climate_Atplus, par_lo = par_lo, par_hi = par_hi,
                             pdf=dnorm, dep_var="MaxGrowth_conv",   initial_temp = 5, temp_red = 0.96, max_iter=100000, hessian = FALSE)

write_results(model5_results_At,"SORTIE parameterization/At3_MaxGrowth_adult.txt", data=F)

#look at residuals
plot(model5_results_At$source_data$predicted, model5_results_At$source_data$MaxGrowth_conv)
SBSmc2_points<-subset(model5_results_At$source_data, model5_results_At$source_data$BGC == "SBSmc2")
points(SBSmc2_points$predicted, SBSmc2_points$MaxGrowth_conv, pch = 2, col = "red")
abline(0,1)

plot(model5_results_At$source_data$aSMR, model5_results_At$source_data$MaxGrowth_conv)
x<-seq(1, 8, .1)
y<-SORTIE.climate.sens(MaxGrowth=model5_results_At$best_pars$MaxGrowth , aL = model5_results_At$best_pars$aL, bL = model5_results_At$best_pars$bL, cL = model5_results_At$best_pars$cL,
                       aH = model5_results_At$best_pars$aH, bH = model5_results_At$best_pars$bH, cH = model5_results_At$best_pars$cH,
                       dL = model5_results_At$best_pars$dL, eL = model5_results_At$best_pars$eL, fL = model5_results_At$best_pars$fL,
                       dH = model5_results_At$best_pars$dH, eH = model5_results_At$best_pars$eH, fH = model5_results_At$best_pars$fH,
                       Xo = model5_results_At$best_pars$Xo, Xb = model5_results_At$best_pars$Xb,
                       SNR = mean(model5_results_At$source_data$SNR), aSMR = x, temp = mean(model5_results_At$source_data$DD5))
lines(x,y)

plot(model5_results_At$source_data$DD5, model5_results_At$source_data$MaxGrowth_conv, xlim = c(300, 1500))
x<-seq(300, 1500, 1)
y<-SORTIE.climate.sens(MaxGrowth=model5_results_At$best_pars$MaxGrowth , aL = model5_results_At$best_pars$aL, bL = model5_results_At$best_pars$bL, cL = model5_results_At$best_pars$cL,
                       aH = model5_results_At$best_pars$aH, bH = model5_results_At$best_pars$bH, cH = model5_results_At$best_pars$cH,
                       dL = model5_results_At$best_pars$dL, eL = model5_results_At$best_pars$eL, fL = model5_results_At$best_pars$fL,
                       dH = model5_results_At$best_pars$dH, eH = model5_results_At$best_pars$eH, fH = model5_results_At$best_pars$fH,
                       Xo = model5_results_At$best_pars$Xo, Xb = model5_results_At$best_pars$Xb,
                       SNR = mean(model5_results_At$source_data$SNR), temp = x, aSMR = mean(model5_results_At$source_data$aSMR))
lines(x,y)


plot(model5_results_At$source_data$SNR, model5_results_At$source_data$MaxGrowth_conv)
x<-seq(1, 5, .1)
y<-SORTIE.climate.sens(MaxGrowth=model5_results_At$best_pars$MaxGrowth , aL = model5_results_At$best_pars$aL, bL = model5_results_At$best_pars$bL, cL = model5_results_At$best_pars$cL,
                       aH = model5_results_At$best_pars$aH, bH = model5_results_At$best_pars$bH, cH = model5_results_At$best_pars$cH,
                       dL = model5_results_At$best_pars$dL, eL = model5_results_At$best_pars$eL, fL = model5_results_At$best_pars$fL,
                       dH = model5_results_At$best_pars$dH, eH = model5_results_At$best_pars$eH, fH = model5_results_At$best_pars$fH,
                       Xo = model5_results_At$best_pars$Xo, Xb = model5_results_At$best_pars$Xb,
                       SNR = x, temp = mean(model5_results_At$source_data$DD5), aSMR = mean(model5_results_At$source_data$aSMR))
lines(x,y)


#MaxGrowth predictions from model. how well do they fit SBSmc2?
SI_conversion$predMaxGrowth <-rep(NA, 4)
SI_conversion$predMaxGrowth[1]<-SORTIE.climate.sens(MaxGrowth=model5_results_At$best_pars$MaxGrowth , aL = model5_results_At$best_pars$aL, bL = model5_results_At$best_pars$bL, cL = model5_results_At$best_pars$cL,
                                                    aH = model5_results_At$best_pars$aH, bH = model5_results_At$best_pars$bH, cH = model5_results_At$best_pars$cH,
                                                    dL = model5_results_At$best_pars$dL, eL = model5_results_At$best_pars$eL, fL = model5_results_At$best_pars$fL,
                                                    dH = model5_results_At$best_pars$dH, eH = model5_results_At$best_pars$eH, fH = model5_results_At$best_pars$fH,
                                                    Xo = model5_results_At$best_pars$Xo, Xb = model5_results_At$best_pars$Xb,
                                                    SNR = 2, temp = 875.73, aSMR = 3)
SI_conversion$predMaxGrowth[2]<-SORTIE.climate.sens(MaxGrowth=model5_results_At$best_pars$MaxGrowth , aL = model5_results_At$best_pars$aL, bL = model5_results_At$best_pars$bL, cL = model5_results_At$best_pars$cL,
                                                    aH = model5_results_At$best_pars$aH, bH = model5_results_At$best_pars$bH, cH = model5_results_At$best_pars$cH,
                                                    dL = model5_results_At$best_pars$dL, eL = model5_results_At$best_pars$eL, fL = model5_results_At$best_pars$fL,
                                                    dH = model5_results_At$best_pars$dH, eH = model5_results_At$best_pars$eH, fH = model5_results_At$best_pars$fH,
                                                    Xo = model5_results_At$best_pars$Xo, Xb = model5_results_At$best_pars$Xb,
                                                    SNR = 3, temp = 875.73, aSMR = 4)
SI_conversion$predMaxGrowth[3]<-SORTIE.climate.sens(MaxGrowth=model5_results_At$best_pars$MaxGrowth , aL = model5_results_At$best_pars$aL, bL = model5_results_At$best_pars$bL, cL = model5_results_At$best_pars$cL,
                                                    aH = model5_results_At$best_pars$aH, bH = model5_results_At$best_pars$bH, cH = model5_results_At$best_pars$cH,
                                                    dL = model5_results_At$best_pars$dL, eL = model5_results_At$best_pars$eL, fL = model5_results_At$best_pars$fL,
                                                    dH = model5_results_At$best_pars$dH, eH = model5_results_At$best_pars$eH, fH = model5_results_At$best_pars$fH,
                                                    Xo = model5_results_At$best_pars$Xo, Xb = model5_results_At$best_pars$Xb,
                                                    SNR = 4, temp = 875.73, aSMR = 5)
SI_conversion$predMaxGrowth[4]<-SORTIE.climate.sens(MaxGrowth=model5_results_At$best_pars$MaxGrowth , aL = model5_results_At$best_pars$aL, bL = model5_results_At$best_pars$bL, cL = model5_results_At$best_pars$cL,
                                                    aH = model5_results_At$best_pars$aH, bH = model5_results_At$best_pars$bH, cH = model5_results_At$best_pars$cH,
                                                    dL = model5_results_At$best_pars$dL, eL = model5_results_At$best_pars$eL, fL = model5_results_At$best_pars$fL,
                                                    dH = model5_results_At$best_pars$dH, eH = model5_results_At$best_pars$eH, fH = model5_results_At$best_pars$fH,
                                                    Xo = model5_results_At$best_pars$Xo, Xb = model5_results_At$best_pars$Xb,
                                                    SNR = 5, temp = 875.73, aSMR = 6)
SI_conversion









###################################################################################################################
###################################################################################################################
####################################         Parameter table           ############################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
















ParameterTable<-matrix(nrow = 4, ncol = 16)
ParameterTable[1,1:16]<-c(model5_results_Sx$R2, model5_results_Sx$best_pars$MaxGrowth, model5_results_Sx$best_pars$aL, model5_results_Sx$best_pars$bL, model5_results_Sx$best_pars$cL,  model5_results_Sx$best_pars$aH, model5_results_Sx$best_pars$bH,
                         model5_results_Sx$best_pars$cH, model5_results_Sx$best_pars$dL, model5_results_Sx$best_pars$eL, model5_results_Sx$best_pars$fL, model5_results_Sx$best_pars$dH, model5_results_Sx$best_pars$eH, model5_results_Sx$best_pars$fH, model5_results_Sx$best_pars$Xo, model5_results_Sx$best_pars$Xb)
ParameterTable[2,1:16]<-c(model5_results_Pl$R2, model5_results_Pl$best_pars$MaxGrowth, model5_results_Pl$best_pars$aL, model5_results_Pl$best_pars$bL, model5_results_Pl$best_pars$cL,  model5_results_Pl$best_pars$aH, model5_results_Pl$best_pars$bH,
                          model5_results_Pl$best_pars$cH, model5_results_Pl$best_pars$dL, model5_results_Pl$best_pars$eL, model5_results_Pl$best_pars$fL, model5_results_Pl$best_pars$dH, model5_results_Pl$best_pars$eH, model5_results_Pl$best_pars$fH, model5_results_Pl$best_pars$Xo, model5_results_Pl$best_pars$Xb)
ParameterTable[3,1:16]<-c(model5_results_Bl$R2, model5_results_Bl$best_pars$MaxGrowth, model5_results_Bl$best_pars$aL, model5_results_Bl$best_pars$bL, model5_results_Bl$best_pars$cL,  model5_results_Bl$best_pars$aH, model5_results_Bl$best_pars$bH,
                          model5_results_Bl$best_pars$cH, model5_results_Bl$best_pars$dL, model5_results_Bl$best_pars$eL, model5_results_Bl$best_pars$fL, model5_results_Bl$best_pars$dH, model5_results_Bl$best_pars$eH, model5_results_Bl$best_pars$fH, model5_results_Bl$best_pars$Xo, model5_results_Bl$best_pars$Xb)
ParameterTable[4,1:16]<-c(model5_results_At$R2, model5_results_At$best_pars$MaxGrowth, model5_results_At$best_pars$aL, model5_results_At$best_pars$bL, model5_results_At$best_pars$cL,  model5_results_At$best_pars$aH, model5_results_At$best_pars$bH,
                          model5_results_At$best_pars$cH, model5_results_At$best_pars$dL, model5_results_At$best_pars$eL, model5_results_At$best_pars$fL, model5_results_At$best_pars$dH, model5_results_At$best_pars$eH, model5_results_At$best_pars$fH, model5_results_At$best_pars$Xo, model5_results_At$best_pars$Xb)
colnames(ParameterTable)<-c( "R2", "MaxGrowth", "aL", "bL", "cL", "aH", "bH", "cH", "dL", "eL", "fL", "dH", "eH", "fH", "Xo", "Xb")
row.names(ParameterTable)<-unique("Sx", "Pl", "Bl", "At")
write.csv(ParameterTable, "ParameterTable_PowerGrowth.csv")
