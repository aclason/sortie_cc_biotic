
#estimating Site index by tree species using aSMR

in_dir <- file.path("01_climate_sensitive_growth","01_data")
out_dir <- file.path("01_climate_sensitive_growth","02_cs_growth_outs")

load(out_dir("rSMR_aSMR_CalcList.RData")) #read in the output of Will's rSMR to aSMR model

Eda <- read.csv(file.path(in_dir,"Edatopic_v10.7.csv"))
sibec <- read.csv(file.path("SIBEC_for_Portfolio2.csv")) ###import actual SI values

colnames(SMRCross) <- c("BGC", "rSMR", "aSMR")
SMRCross$rSMR <- gsub("[[:alpha:]]","", SMRCross$rSMR)
SMRCross$aSMR <- round(SMRCross$aSMR, digits = 0)

sibec <- subset(sibec, sibec$PlotCountSpp>0 )## remove site series with no plots

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

##############below code has already been used to prep data for each species.
#It needs to be re-examined in future iterations
sibec <- sibec[sibec$TreeSpp == "Pl",]##choose species
sibec <- sibec[!is.na(sibec$BGCUnit),]

#Erica note: as is, this is giving the number of site series listed in "sibec" for each BGCUnit
#   I don't think this is what we want because MeanPlotSiteIndex also has si values with no plots
#   listed, which must be guesses?
#numPlots <- aggregate(MeanPlotSiteIndex ~ BGCUnit, sibec, FUN = length)

###remove SI values in zones with little data (use 2 for Lw and Py)
#BGC <- numPlots$BGCUnit[numPlots$MeanPlotSiteIndex > 4]
#sibec <- sibec[sibec$BGCUnit %in% BGC,]

#instead, lets first remove site series that don't have any plots...this is done above in the code
sibec_plots<-subset(sibec, sibec$PlotCountSpp > 0 )
#remaining site series have at least 7 plots with data (for Bl), that seems sufficient to include
min(sibec_plots$PlotCountSpp)
# then we can use above code but change numPlots to numSiteSeries which describes the
# object more accurately
numSiteSeries <- aggregate(SiteSeries ~ BGCUnit, sibec_plots, FUN = length)
#changed to >= 4 for Bl so enough BGCs would be included
BGC <- numSiteSeries$BGCUnit[numSiteSeries$SiteSeries >= 4]
sibec <- sibec_plots[sibec_plots$BGCUnit %in% BGC,]

####for BL, 3 BGC  unit site series don't match up between sibec data and edatope so
# we are losing data during this merge need another cross walk table here for "ESSFdc1"
# "ESSFwc1" and "SBSwk3"(doesn't have data for edatopic for SBSwk3/08 and another site series)
# SBSwk3 doesn't get through the merge while looped but it does individually because the Edatopic
# is missing for 2 site series and there is a is.na function

#try converting sibec$BGCUnit to character to fix loop error
sibec$BGCUnit <- as.character(sibec$BGCUnit)
Eda$MergedBGC <- as.character(Eda$MergedBGC)
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
