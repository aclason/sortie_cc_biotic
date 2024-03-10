rm(list=ls())   #cleans the workspace so all previous objects are deleted

#This script brings in PSP data and takes samples of interest for SORTIE-climate paper
#then creates multi-variate dataset with diameter distribution to group stands for sample selection

setwd("C:/Users/elilles/OneDrive - Government of BC/SORTIE_climate_scenarios")
#load libraries
library(readxl)
library(plyr)
library(ecodist)

  ##############################################
#######      bring in PSP data for SBSmc2 ################
   ##############################################
setwd("C:/Users/elilles/OneDrive - Government of BC/SORTIE ESSF parameterization/allometry")

#Include SBSmc subzones and variants
IncludeBEC<-c("SBSmc2", "SBSmc1", "SBSmc3")
#first bring in plot level data to identify which plots to sample from
#TSA_03 (Bulkley)
TSA_03_sample<-read.csv("datasets/PSPdata_TSA03_sample.csv")
str(TSA_03_sample)
TSA_03_sample<-subset(TSA_03_sample, TSA_03_sample$beclabel %in% IncludeBEC)
TSA_03_SAMP_ID<-unique(TSA_03_sample$SAMP_ID)

#TSA_14 (Lakes)
TSA_14_sample<-read_excel(path = "datasets/TSA14.xlsx", sheet = "sample")
TSA_14_sample<-as.data.frame(TSA_14_sample)
TSA_14_sample<-subset(TSA_14_sample, TSA_14_sample$beclabel %in% IncludeBEC)
TSA_14_SAMP_ID<-unique(TSA_14_sample$SAMP_ID)

#TSA_20 (Morice)
TSA_20_sample<-read_excel(path = "datasets/TSA20.xlsx", sheet = "sample")
TSA_20_sample<-as.data.frame(TSA_20_sample)
TSA_20_sample<-subset(TSA_20_sample, TSA_20_sample$beclabel %in% IncludeBEC)
TSA_20_SAMP_ID<-unique(TSA_20_sample$SAMP_ID)

#TSA_24 (PrinceGeorge)
TSA_24_sample<-read.csv("datasets/PSPdata_TSA24_sample.csv")
str(TSA_24_sample)
TSA_24_sample<-subset(TSA_24_sample, TSA_24_sample$beclabel %in% IncludeBEC)
TSA_24_SAMP_ID<-unique(TSA_24_sample$SAMP_ID)

#TSA_26 (Quesnel)
TSA_26_sample<-read_excel(path = "datasets/TSA26.xlsx", sheet = "sample")
TSA_26_sample<-as.data.frame(TSA_26_sample)
TSA_26_sample<-subset(TSA_26_sample, TSA_26_sample$beclabel %in% IncludeBEC)
TSA_26_SAMP_ID<-unique(TSA_26_sample$SAMP_ID)

#then bring in tree data and select trees from the right plots 
TSA_03_tree<-read.csv("datasets/PSPdata_TSA03_tree.csv")
str(TSA_03_tree)
TSA_03_tree<-subset(TSA_03_tree, TSA_03_tree$samp_id %in% TSA_03_SAMP_ID)

TSA_14_tree<-read_excel(path = "datasets/TSA14.xlsx", sheet = "tree")
TSA_14_tree<-as.data.frame(TSA_14_tree)
TSA_14_tree<-subset(TSA_14_tree, TSA_14_tree$samp_id %in% TSA_14_SAMP_ID)

TSA_20_tree<-read_excel(path = "datasets/TSA20.xlsx", sheet = "tree")
TSA_20_tree<-as.data.frame(TSA_20_tree)
TSA_20_tree<-subset(TSA_20_tree, TSA_20_tree$samp_id %in% TSA_20_SAMP_ID)

TSA_24_tree<-read.csv("datasets/PSPdata_TSA24_tree.csv")
str(TSA_24_tree)
TSA_24_tree<-subset(TSA_24_tree, TSA_24_tree$samp_id %in% TSA_24_SAMP_ID)

TSA_26_tree<-read_excel(path = "datasets/TSA26.xlsx", sheet = "tree")
TSA_26_tree<-as.data.frame(TSA_26_tree)
TSA_26_tree<-subset(TSA_26_tree, TSA_26_tree$samp_id %in% TSA_26_SAMP_ID)

#bind three TSA datasets together, use sample with a 4 cm dbh compile limit.. available for all plots, only G and t have 2 cm dbh limit
SBSmc_sample<- rbind(TSA_03_sample, TSA_14_sample, TSA_20_sample, TSA_26_sample, TSA_24_sample)
#use the most recent measurement only
SBSmc_sample<-subset(SBSmc_sample, SBSmc_sample$meas_last == "Y")

#plots in samples
table(SBSmc_sample$samp_id, SBSmc_sample$dbhlimit_compile)
table(SBSmc_sample$sampletype, SBSmc_sample$dbhlimit_compile)

SBSmc_sample<-subset(SBSmc_sample, SBSmc_sample$dbhlimit_compile == 4)
SBSmc_tree<- rbind(TSA_03_tree, TSA_14_tree, TSA_20_tree, TSA_24_tree, TSA_26_tree)
str(SBSmc_tree)

#fixed and variable radius plots included in data
#VRI data is all variable radius except on plot. will remove that plot
table(SBSmc_sample$sampletype, SBSmc_sample$plot_typ)
subset(SBSmc_sample, SBSmc_sample$sampletype == "VRI" & SBSmc_sample$plot_typ == "F")
SBSmc_sample<-subset(SBSmc_sample, SBSmc_sample$SAMP_ID != "0031_0113_VRI")

#clean data
#judgement call here:
#only use the last measurement of each plot so that we are projecting forward recent stands 
#and have the best data, including the smallest size classes
#consider limiting sample to only those with small tree counts
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$meas_last == "Y")  

unique(SBSmc_tree$ld)
SBSmc_tree[which(SBSmc_tree$ld == "C"),] #C is for cut trees, not sure what to do because there in only one I'll leave it in a live
live_list<-c( "L",  "I" , "V" , "C") 
SBSmc_tree$LiveDead<-ifelse(SBSmc_tree$ld %in% live_list, "L", "D")

#I don't know why there are trees without DBHs, they look like they were used for site index so maybe they 
#were outside of the plot, judgement call to remove them
SBSmc_tree[is.na(SBSmc_tree$dbh) == TRUE,] 
SBSmc_tree<- SBSmc_tree[is.na(SBSmc_tree$dbh) == FALSE,] 
hist(SBSmc_tree$dbh)
min(SBSmc_tree$dbh)

#remove trees below compile limit
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$dbh >=4)


unique(SBSmc_tree$species)
#where are weird species coming from?

SBSmc_tree[which(SBSmc_tree$species == "DR"),]
#remove plots with red alder (DR)
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "0031_0097_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "66048 G000177")

SBSmc_tree[which(SBSmc_tree$species == "CW"),] #only one tree, OK to leave in
SBSmc_tree[which(SBSmc_tree$species == "HW"),] 
#remove plots with more than one hemlock (HW)
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "0031_0122_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "0031_0125_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "67011 G000054")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "67022 G000052")

SBSmc_tree[which(SBSmc_tree$species == "DM"),] #mountain alder ... ok to keep in, can lump with deciduous
SBSmc_tree[which(SBSmc_tree$species == "FD"),] #only two FD, can keep plots in

SBSmc_tree[which(SBSmc_tree$species == "SB"),]
#remove plots with more than one black spruce, not the type of stand we are trying to model
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "0031_0108_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "0031_0130_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "0141_0076_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "0141_0077_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "0141_0102_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "66004 G000027")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "66014 G000120")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "024C_5596_CMI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "CMI4_0298_NFI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "DPG1_0028_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "DVA1_0021_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "DVA1_0023_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "DVA1_0029_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "DVA1_0049_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "DVA1_0060_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "DVA1_0061_VRI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "DVA1_0066_VRI")

SBSmc_tree[which(SBSmc_tree$species == "BA"),] 
#remove plots with more than one BA, 
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "024C_7416_CMI")

SBSmc_tree[which(SBSmc_tree$species == "W"),] #Willow, can lump with deciduous
#some plot seems to be empty... will remove
SBSmc_tree[which(SBSmc_tree$species == ""),] 
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "0031_0105_VRI")
SBSmc_tree[is.na(SBSmc_tree$species) == TRUE,] 
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "0142_1171_CMI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "0142_8696_CMI")
SBSmc_tree<-subset(SBSmc_tree, SBSmc_tree$samp_id != "BFP1_0025_VRI")

#lumping species codes
SBSmc_tree$species[which(SBSmc_tree$species == "B")]<- "BL" #assuming these are BL
SBSmc_tree$species[which(SBSmc_tree$species == "S")]<- "SX" #assuming these are SX
SBSmc_tree$species[which(SBSmc_tree$species == "SW")]<- "SX" #assuming these are SX

unique(SBSmc_tree$species)
deciduous<-c("AC", "AT", "EP", "DM", "W")
SBSconifer<-c("SX", "BL", "PL")
SBSmc_tree$SppGrp<- ifelse(SBSmc_tree$species == "SX", "SX",
                     ifelse(SBSmc_tree$species == "PL", "PL",
                     ifelse( SBSmc_tree$species == "BL", "BL",
                    ifelse( SBSmc_tree$species %in% deciduous, "DECID", "OTHER"))))
SBSmc_tree$SppGrp<-ifelse(SBSmc_tree$species %in%  SBSconifer & SBSmc_tree$LiveDead == "D", "CON_D",SBSmc_tree$SppGrp)
SBSmc_tree$SppGrp<-ifelse(SBSmc_tree$species %in%  deciduous & SBSmc_tree$LiveDead == "D", "DECID_D",SBSmc_tree$SppGrp)

 ######################
#look into when these plots were measured and only one measurement per plot
table(SBSmc_tree$samp_id, SBSmc_tree$meas_yr)
############################################################
 #splitting data into DBH classes

##########
###########Functions
###########

#loop to create size classes to get data into
#dropping size class of 5 because of 4 cm dbh limit for CMI, VRI and YSM plots
sizeClasses <- c(  10, 20, 30, 80)
inits <- vector()
for(i in 1:length(sizeClasses)){
  inits[i] <- paste0("Init.Dens.",sizeClasses[i])
}
init.values <- matrix(nrow=length(sizeClasses),ncol=2)
row.names(init.values) <- inits
for(j in 1:length(sizeClasses)){
  init.values[j,1] <- sizeClasses[j]-5
  init.values[j,2] <- sizeClasses[j]
}
init.values
init.values[1,1]<-4
init.values[2,1]<-10
init.values[3,1]<-20
init.values[4,1]<-30
init.values

###in this code the lower limit is included in the group and the upper limit is in next the group above
dbhClassphf <- function(dbhSPH){ ###calulates SPH in DBH categories
  ###cutoffs for size classes are stored in init.values
  results <- vector("numeric",2)
  for(i in 1:length(sizeClasses)){
    classSum <- sum(dbhSPH$phf_tree[dbhSPH$dbh >= init.values[i,1] & dbhSPH$dbh < init.values[i,2]])
    results[i] <- classSum
  }
  return(results)
}
#########
#######Stems per ha of different species and live dead categories
######

#########
#######prep all data together but remove VRI which has multiple plots, then process VRI separately and bind together again
######

length(unique(SBSmc_tree$samp_id))
labels<-ddply(SBSmc_tree[c("samp_id",  "phf_tree")], .(samp_id), numcolwise(sum))
labels<-labels[c("samp_id")]
PSPs_SppGrp<-unique(SBSmc_tree$SppGrp)

DECIDtemp <-subset(SBSmc_tree, SBSmc_tree$SppGrp == "DECID" )
SBSmc_tree_output_DECID<-ddply(DECIDtemp[c("samp_id", "dbh", "phf_tree")], .(samp_id), dbhClassphf) #this is like a pivot table in excel
colnames(SBSmc_tree_output_DECID)<-c("samp_id",  inits)
#get labels to make sure plots with zero values are entered
SBSmc_tree_output_DECID<-merge(labels, SBSmc_tree_output_DECID, by = c("samp_id"), all.x = TRUE)
SBSmc_tree_output_DECID[is.na(SBSmc_tree_output_DECID)] <- 0
SBSmc_tree_output_DECID$DECID_sm<- SBSmc_tree_output_DECID$Init.Dens.10
SBSmc_tree_output_DECID$DECID_lg<-SBSmc_tree_output_DECID$Init.Dens.20 + SBSmc_tree_output_DECID$Init.Dens.30 + SBSmc_tree_output_DECID$Init.Dens.80
SBSmc_tree_output_DECID<-SBSmc_tree_output_DECID[c("samp_id", "DECID_sm", "DECID_lg")]

PLtemp <-subset(SBSmc_tree, SBSmc_tree$SppGrp == "PL" )
SBSmc_tree_output_PL<-ddply(PLtemp[c("samp_id","dbh", "phf_tree")], .(samp_id), dbhClassphf) #this is like a pivot table in excel
colnames(SBSmc_tree_output_PL)<-c("samp_id", "PLInit.Dens.10", "PLInit.Dens.20", "PLInit.Dens.30", "PLInit.Dens.80")
#get labels to make sure plots with zero values are entered
SBSmc_tree_output_PL<-merge(labels, SBSmc_tree_output_PL, by = c("samp_id"), all.x = TRUE)
SBSmc_tree_output_PL[is.na(SBSmc_tree_output_PL)] <- 0

BLtemp <-subset(SBSmc_tree, SBSmc_tree$SppGrp == "BL" )
SBSmc_tree_output_BL<-ddply(BLtemp[c("samp_id","dbh", "phf_tree")], .(samp_id), dbhClassphf) #this is like a pivot table in excel
colnames(SBSmc_tree_output_BL)<-c("samp_id", "BLInit.Dens.10", "BLInit.Dens.20", "BLInit.Dens.30", "BLInit.Dens.80")
#get labels to make sure plots with zero values are entered
SBSmc_tree_output_BL<-merge(labels, SBSmc_tree_output_BL, by =  c("samp_id"), all.x = TRUE)
SBSmc_tree_output_BL[is.na(SBSmc_tree_output_BL)] <- 0

SXtemp <-subset(SBSmc_tree, SBSmc_tree$SppGrp == "SX" )
SBSmc_tree_output_SX<-ddply(SXtemp[c("samp_id","dbh", "phf_tree")], .(samp_id), dbhClassphf) #this is like a pivot taSXe in excel
colnames(SBSmc_tree_output_SX)<-c("samp_id", "SXInit.Dens.10", "SXInit.Dens.20", "SXInit.Dens.30", "SXInit.Dens.80")
#get labels to make sure plots with zero values are entered
SBSmc_tree_output_SX<-merge(labels, SBSmc_tree_output_SX, by =  c("samp_id"), all.x = TRUE)
SBSmc_tree_output_SX[is.na(SBSmc_tree_output_SX)] <- 0

DECID_Dtemp <-subset(SBSmc_tree, SBSmc_tree$SppGrp == "DECID_D" )
SBSmc_tree_output_DECID_D<-ddply(DECID_Dtemp[c("samp_id","dbh", "phf_tree")], .(samp_id), dbhClassphf) #this is like a pivot taDECID_De in excel
colnames(SBSmc_tree_output_DECID_D)<-c("samp_id", "DECID_DInit.Dens.10", "DECID_DInit.Dens.20", "DECID_DInit.Dens.30", "DECID_DInit.Dens.80")
#get labels to make sure plots with zero values are entered
SBSmc_tree_output_DECID_D<-merge(labels, SBSmc_tree_output_DECID_D, by =  c("samp_id"), all.x = TRUE)
SBSmc_tree_output_DECID_D[is.na(SBSmc_tree_output_DECID_D)] <- 0
SBSmc_tree_output_DECID_D$DECID_D<-SBSmc_tree_output_DECID_D$DECID_DInit.Dens.10+
  SBSmc_tree_output_DECID_D$DECID_DInit.Dens.20 +SBSmc_tree_output_DECID_D$DECID_DInit.Dens.30 +SBSmc_tree_output_DECID_D$DECID_DInit.Dens.80
SBSmc_tree_output_DECID_D<- SBSmc_tree_output_DECID_D[c("samp_id", "DECID_D")]

CON_Dtemp <-subset(SBSmc_tree, SBSmc_tree$SppGrp == "CON_D" )
SBSmc_tree_output_CON_D<-ddply(CON_Dtemp[c("samp_id","dbh", "phf_tree")], .(samp_id), dbhClassphf) #this is like a pivot table in excel
colnames(SBSmc_tree_output_CON_D)<-c("samp_id", "CON_DInit.Dens.10", "CON_DInit.Dens.20", "CON_DInit.Dens.30", "CON_DInit.Dens.80")
#get labels to make sure plots with zero values are entered
SBSmc_tree_output_CON_D<-merge(labels, SBSmc_tree_output_CON_D, by =  c("samp_id"), all.x = TRUE)
SBSmc_tree_output_CON_D[is.na(SBSmc_tree_output_CON_D)] <- 0
SBSmc_tree_output_CON_D$CON_D_sm<-SBSmc_tree_output_CON_D$CON_DInit.Dens.10
SBSmc_tree_output_CON_D$CON_D_lg<-SBSmc_tree_output_CON_D$CON_DInit.Dens.20 + SBSmc_tree_output_CON_D$CON_DInit.Dens.30 + SBSmc_tree_output_CON_D$CON_DInit.Dens.80
SBSmc_tree_output_CON_D<-SBSmc_tree_output_CON_D[c("samp_id", "CON_D_sm", "CON_D_lg")]


SBSmc_tree_output<-merge(SBSmc_tree_output_PL, SBSmc_tree_output_SX, by =  c("samp_id"))
SBSmc_tree_output<-merge(SBSmc_tree_output, SBSmc_tree_output_BL, by =  c("samp_id"))
SBSmc_tree_output<-merge(SBSmc_tree_output, SBSmc_tree_output_DECID, by =  c("samp_id"))
SBSmc_tree_output<-merge(SBSmc_tree_output, SBSmc_tree_output_DECID_D, by =  c("samp_id"))
SBSmc_tree_output<-merge(SBSmc_tree_output, SBSmc_tree_output_CON_D, by =  c("samp_id"))
head(SBSmc_tree_output)
#remove VRI samples
SBSmc_tree_output.noVRI<-subset(SBSmc_tree_output, substr(SBSmc_tree_output$samp_id, 11,13) != "VRI")
unique(SBSmc_tree_output.noVRI$samp_id)

#checking, checks out without VRI
SBSmc_tree_output.noVRI$SPHsums<-rowSums(SBSmc_tree_output.noVRI[,2:length(SBSmc_tree_output.noVRI)])
SBSmc_sample$samp_id<- SBSmc_sample$SAMP_ID
test<-merge(SBSmc_tree_output.noVRI, SBSmc_sample[c("samp_id", "sampletype", "stemsha_LIV", "stemsha_DP", "stemsha_DU", "stemsha_DF")], by = "samp_id")
test$totalSPH<-test$stemsha_LIV+ test$stemsha_DP+ test$stemsha_DU +test$stemsha_DF
plot(test$totalSPH, test$SPHsums)
abline(0,1)


#########
#######prep VRI which has multiple plots, then process VRI separately and bind together again
######
SBSmc_tree.VRI<-subset(SBSmc_tree, substr(SBSmc_tree$samp_id, 11,13) == "VRI")
length(unique(SBSmc_tree.VRI$samp_id))
labels<-ddply(SBSmc_tree.VRI[c("samp_id", "plot_no", "phf_tree")], .(samp_id, plot_no), numcolwise(sum))
labels<-labels[c("samp_id", "plot_no")]
PSPs_SppGrp<-unique(SBSmc_tree.VRI$SppGrp)

DECIDtemp <-subset(SBSmc_tree.VRI, SBSmc_tree.VRI$SppGrp == "DECID" )
SBSmc_tree.VRI_output_DECID<-ddply(DECIDtemp[c("samp_id", "plot_no","dbh", "phf_tree")], .(samp_id, plot_no), dbhClassphf) #this is like a pivot table in excel
colnames(SBSmc_tree.VRI_output_DECID)<-c("samp_id", "plot_no", inits)
#get labels to make sure plots with zero values are entered
SBSmc_tree.VRI_output_DECID<-merge(labels, SBSmc_tree.VRI_output_DECID, by = c("samp_id", "plot_no"), all.x = TRUE)
SBSmc_tree.VRI_output_DECID[is.na(SBSmc_tree.VRI_output_DECID)] <- 0
SBSmc_tree.VRI_output_DECID$DECID_sm<- SBSmc_tree.VRI_output_DECID$Init.Dens.10
SBSmc_tree.VRI_output_DECID$DECID_lg<-SBSmc_tree.VRI_output_DECID$Init.Dens.20 + SBSmc_tree.VRI_output_DECID$Init.Dens.30 + SBSmc_tree.VRI_output_DECID$Init.Dens.80
SBSmc_tree.VRI_output_DECID<-SBSmc_tree.VRI_output_DECID[c("samp_id","plot_no","DECID_sm", "DECID_lg")]

PLtemp <-subset(SBSmc_tree.VRI, SBSmc_tree.VRI$SppGrp == "PL" )
SBSmc_tree.VRI_output_PL<-ddply(PLtemp[c("samp_id","plot_no","dbh", "phf_tree")], .(samp_id, plot_no), dbhClassphf) #this is like a pivot table in excel
colnames(SBSmc_tree.VRI_output_PL)<-c("samp_id", "plot_no", "PLInit.Dens.10", "PLInit.Dens.20", "PLInit.Dens.30", "PLInit.Dens.80")
#get labels to make sure plots with zero values are entered
SBSmc_tree.VRI_output_PL<-merge(labels, SBSmc_tree.VRI_output_PL, by = c("samp_id", "plot_no"), all.x = TRUE)
SBSmc_tree.VRI_output_PL[is.na(SBSmc_tree.VRI_output_PL)] <- 0

BLtemp <-subset(SBSmc_tree.VRI, SBSmc_tree.VRI$SppGrp == "BL" )
SBSmc_tree.VRI_output_BL<-ddply(BLtemp[c("samp_id","plot_no","dbh", "phf_tree")], .(samp_id, plot_no), dbhClassphf) #this is like a pivot table in excel
colnames(SBSmc_tree.VRI_output_BL)<-c("samp_id","plot_no", "BLInit.Dens.10", "BLInit.Dens.20", "BLInit.Dens.30", "BLInit.Dens.80")
#get labels to make sure plots with zero values are entered
SBSmc_tree.VRI_output_BL<-merge(labels, SBSmc_tree.VRI_output_BL, by =  c("samp_id", "plot_no"), all.x = TRUE)
SBSmc_tree.VRI_output_BL[is.na(SBSmc_tree.VRI_output_BL)] <- 0

SXtemp <-subset(SBSmc_tree.VRI, SBSmc_tree.VRI$SppGrp == "SX" )
SBSmc_tree.VRI_output_SX<-ddply(SXtemp[c("samp_id","plot_no","dbh", "phf_tree")], .(samp_id, plot_no), dbhClassphf) #this is like a pivot taSXe in excel
colnames(SBSmc_tree.VRI_output_SX)<-c("samp_id","plot_no", "SXInit.Dens.10", "SXInit.Dens.20", "SXInit.Dens.30", "SXInit.Dens.80")
#get labels to make sure plots with zero values are entered
SBSmc_tree.VRI_output_SX<-merge(labels, SBSmc_tree.VRI_output_SX, by =  c("samp_id", "plot_no"), all.x = TRUE)
SBSmc_tree.VRI_output_SX[is.na(SBSmc_tree.VRI_output_SX)] <- 0

DECID_Dtemp <-subset(SBSmc_tree.VRI, SBSmc_tree.VRI$SppGrp == "DECID_D" )
SBSmc_tree.VRI_output_DECID_D<-ddply(DECID_Dtemp[c("samp_id","plot_no","dbh", "phf_tree")], .(samp_id, plot_no), dbhClassphf) #this is like a pivot taDECID_De in excel
colnames(SBSmc_tree.VRI_output_DECID_D)<-c("samp_id","plot_no", "DECID_DInit.Dens.10", "DECID_DInit.Dens.20", "DECID_DInit.Dens.30", "DECID_DInit.Dens.80")
#get labels to make sure plots with zero values are entered
SBSmc_tree.VRI_output_DECID_D<-merge(labels, SBSmc_tree.VRI_output_DECID_D, by =  c("samp_id", "plot_no"), all.x = TRUE)
SBSmc_tree.VRI_output_DECID_D[is.na(SBSmc_tree.VRI_output_DECID_D)] <- 0
SBSmc_tree.VRI_output_DECID_D$DECID_D<-SBSmc_tree.VRI_output_DECID_D$DECID_DInit.Dens.10+
  SBSmc_tree.VRI_output_DECID_D$DECID_DInit.Dens.20 +SBSmc_tree.VRI_output_DECID_D$DECID_DInit.Dens.30 +SBSmc_tree.VRI_output_DECID_D$DECID_DInit.Dens.80
SBSmc_tree.VRI_output_DECID_D<- SBSmc_tree.VRI_output_DECID_D[c("samp_id","plot_no", "DECID_D")]

CON_Dtemp <-subset(SBSmc_tree.VRI, SBSmc_tree.VRI$SppGrp == "CON_D" )
SBSmc_tree.VRI_output_CON_D<-ddply(CON_Dtemp[c("samp_id","plot_no","dbh", "phf_tree")], .(samp_id, plot_no), dbhClassphf) #this is like a pivot table in excel
colnames(SBSmc_tree.VRI_output_CON_D)<-c("samp_id","plot_no", "CON_DInit.Dens.10", "CON_DInit.Dens.20", "CON_DInit.Dens.30", "CON_DInit.Dens.80")
#get labels to make sure plots with zero values are entered
SBSmc_tree.VRI_output_CON_D<-merge(labels, SBSmc_tree.VRI_output_CON_D, by =  c("samp_id", "plot_no"), all.x = TRUE)
SBSmc_tree.VRI_output_CON_D[is.na(SBSmc_tree.VRI_output_CON_D)] <- 0
SBSmc_tree.VRI_output_CON_D$CON_D_sm<- SBSmc_tree.VRI_output_CON_D$CON_DInit.Dens.10
SBSmc_tree.VRI_output_CON_D$CON_D_lg<-SBSmc_tree.VRI_output_CON_D$CON_DInit.Dens.20 + SBSmc_tree.VRI_output_CON_D$CON_DInit.Dens.30 + SBSmc_tree.VRI_output_CON_D$CON_DInit.Dens.80
SBSmc_tree.VRI_output_CON_D<-SBSmc_tree.VRI_output_CON_D[c("samp_id","plot_no", "CON_D_sm", "CON_D_lg")]


SBSmc_tree.VRI_output<-merge(SBSmc_tree.VRI_output_PL, SBSmc_tree.VRI_output_SX, by =  c("samp_id", "plot_no"))
SBSmc_tree.VRI_output<-merge(SBSmc_tree.VRI_output, SBSmc_tree.VRI_output_BL, by =  c("samp_id", "plot_no"))
SBSmc_tree.VRI_output<-merge(SBSmc_tree.VRI_output, SBSmc_tree.VRI_output_DECID, by =  c("samp_id", "plot_no"))
SBSmc_tree.VRI_output<-merge(SBSmc_tree.VRI_output, SBSmc_tree.VRI_output_DECID_D, by =  c("samp_id", "plot_no"))
SBSmc_tree.VRI_output<-merge(SBSmc_tree.VRI_output, SBSmc_tree.VRI_output_CON_D, by =  c("samp_id", "plot_no"))
head(SBSmc_tree.VRI_output)
SBSmc_tree.VRI_output_s<-ddply(SBSmc_tree.VRI_output, .(samp_id), numcolwise(mean)) #this is like a pivot table in excel

#checking, checks out (more or less, possible errors in data)
SBSmc_tree.VRI_output_s$SPHsums<-rowSums(SBSmc_tree.VRI_output_s[,2:length(SBSmc_tree.VRI_output_s)])
test2<-merge(SBSmc_tree.VRI_output_s, SBSmc_sample[c("samp_id", "sampletype", "stemsha_LIV", "stemsha_DP", "stemsha_DU", "stemsha_DF")], by = "samp_id")
test2$totalSPH<-test2$stemsha_LIV+ test2$stemsha_DP+ test2$stemsha_DU +test2$stemsha_DF
plot(test2$totalSPH, test2$SPHsums)
abline(0,1)

#bind together
SBSmc_tree_output<-rbind(SBSmc_tree_output.noVRI, SBSmc_tree.VRI_output_s)
SBSmc_tree_output$SPHsums <-NULL #after checking, no need for this column
#########
####### Running NMDS
######

#make the rownames the Tree.ID
row.names(SBSmc_tree_output)<-SBSmc_tree_output$samp_id
Codes<-SBSmc_tree_output[,1] #save the info for each tree in a codes file

#get bray-curtis distances, run NMDS and output minimum distance configuration with stress and R2 info
dist<-distance(SBSmc_tree_output[,2:length(SBSmc_tree_output)], "bray-curtis")
distNMDS<-nmds(dist, mindim = 2, maxdim = 2, maxit = 100000)
dist.scors<-nmds.min(distNMDS)

#create dataframe with nmds results and attach to codes
NMS_labelled_points<- data.frame(dist.scors, Codes )

setwd("C:/Users/elilles/OneDrive - Government of BC/SORTIE_climate_scenarios")
#write.csv(NMS_labelled_points, file = "NMS_results_PSP_SBSmc_all_1.csv")

#plot nmds points with symbols for different categories
windows()
NMS_labelled_points$samp_id<-NMS_labelled_points$Codes


#drop 16 plots with no age
test<-SBSmc_sample[is.na(SBSmc_sample$tot_stand_age),]
length(unique(test$SAMP_ID))
SBSmc_sample<- SBSmc_sample[is.na(SBSmc_sample$tot_stand_age) == FALSE,]

#drop 13 plots with no site index or negative site index
test.SI<-SBSmc_sample[is.na(SBSmc_sample$lead_si1),]
length(unique(test.SI$SAMP_ID))
SBSmc_sample<- subset( SBSmc_sample, SBSmc_sample$lead_si1 >0)
                       
#bring in info from sample data... don't keep all of x or y so that we lose samples that were excluded for one reason or another
NMS_labelled_points_l<-merge(NMS_labelled_points, SBSmc_sample[c("samp_id", "sampletype","tsa", "meas_yr",  "spc_live1", "tot_stand_age", "dbhq_LIV", "baha_LIV", "stemsha_LIV", "wsvha_LIV", "lead_si1")], by = "samp_id")
NMS_labelled_points_l<-merge(NMS_labelled_points_l, SBSmc_tree_output, by = "samp_id")
Age1<-subset(NMS_labelled_points_l, NMS_labelled_points_l$tot_stand_age<= 20)
Age2<-subset(NMS_labelled_points_l, NMS_labelled_points_l$tot_stand_age> 20 &NMS_labelled_points_l$tot_stand_age<= 40)
young_mat<-subset(NMS_labelled_points_l, NMS_labelled_points_l$tot_stand_age> 40&NMS_labelled_points_l$tot_stand_age<= 80)
Age5<-subset(NMS_labelled_points_l, NMS_labelled_points_l$tot_stand_age> 80&NMS_labelled_points_l$tot_stand_age<= 100)
mature<-subset(NMS_labelled_points_l, NMS_labelled_points_l$tot_stand_age> 100&NMS_labelled_points_l$tot_stand_age<= 140)
old<-subset(NMS_labelled_points_l, NMS_labelled_points_l$tot_stand_age> 140&NMS_labelled_points_l$tot_stand_age<= 230)
v.old<-subset(NMS_labelled_points_l, NMS_labelled_points_l$tot_stand_age> 230&NMS_labelled_points_l$tot_stand_age< 999)

plot(x=NMS_labelled_points_l$X1, y=NMS_labelled_points_l$X2, xlim = c(-.7,.7), ylim = c(-.7,.7), xaxt = "n", yaxt = "n", xlab = "", ylab = "", pch = 16)
points(x=Age2$X1, y=Age2$X2,pch = 10, col = "red")
points(x=young_mat$X1, y=young_mat$X2,pch = 10,  col = "orange")
points(x=Age5$X1, y=Age5$X2,pch = 2,  col = "green")
points(x=old$X1, y=old$X2,pch = 3, col = "blue")
points(x=v.old$X1, y=v.old$X2,pch = 5, col = "purple")
legend(x="bottomright", legend = c( "Age2", "40-80", 'Age5','140-230','>230'), pch = c(10,10,2,3,5), col = c("red", "orange", "green", "blue", "purple"))

#add vectors
vectors=vf(dist.scors, SBSmc_tree_output[,2:length(SBSmc_tree_output)], nperm=10)
#vectors<-subset(vectors, vectors[,4] <= 0.1) trying to only plot vectors with p-values of 0.1 or lower
plot(vectors, 
     #len=0.1,
     col="red")
####

#show which plots are circum SBSmc2 mesic (SBSmc2 01 site index is 17.9 for Spruce and 18.8 for Pine)

circum_mesic<-subset(NMS_labelled_points_l, NMS_labelled_points_l$lead_si1 >17 & NMS_labelled_points_l$lead_si1<= 19.5 )
richer<-subset(NMS_labelled_points_l, NMS_labelled_points_l$lead_si1 > 19.5 )
poorer<-subset(NMS_labelled_points_l, NMS_labelled_points_l$lead_si1 <= 17 )

plot(x=NMS_labelled_points_l$X1, y=NMS_labelled_points_l$X2, xlim = c(-.7,.7), ylim = c(-.7,.7), xaxt = "n", yaxt = "n", xlab = "", ylab = "", pch = 16)
points(x=circum_mesic$X1, y=circum_mesic$X2,pch = 10, col = "red")
points(x=richer$X1, y=richer$X2,pch = 10,  col = "orange")
points(x=poorer$X1, y=poorer$X2,pch = 5,  col = "green")
legend(x="bottomright", legend = c( "circummesic", "richer", 'poorer'), pch = c(10,10,2), col = c("red", "orange", "green"))

##### Leading species
NMS_labelled_points_l$PinePer <-100*rowSums(NMS_labelled_points_l[c("PLInit.Dens.10", "PLInit.Dens.20", "PLInit.Dens.30","PLInit.Dens.80")])/  NMS_labelled_points_l$stemsha_LIV
NMS_labelled_points_l$SprucePer <-100*rowSums(NMS_labelled_points_l[c("SXInit.Dens.10", "SXInit.Dens.20", "SXInit.Dens.30","SXInit.Dens.80")])/  NMS_labelled_points_l$stemsha_LIV
NMS_labelled_points_l$FirPer <-100*rowSums(NMS_labelled_points_l[c("BLInit.Dens.10", "BLInit.Dens.20", "BLInit.Dens.30","BLInit.Dens.80")])/  NMS_labelled_points_l$stemsha_LIV
NMS_labelled_points_l$DecidPer <-100*rowSums(NMS_labelled_points_l[c("DECID_sm",  "DECID_lg" )])/  NMS_labelled_points_l$stemsha_LIV

PL_dominant<-subset(NMS_labelled_points_l, NMS_labelled_points_l$PinePer >= 70 & NMS_labelled_points_l$PinePer <=100)
SX_dominant<-subset(NMS_labelled_points_l, NMS_labelled_points_l$SprucePer >= 70 & NMS_labelled_points_l$SprucePer <=100 )
BL_dominant<-subset(NMS_labelled_points_l, NMS_labelled_points_l$FirPer >= 70 & NMS_labelled_points_l$FirPer <=100)
DECID_mix<-subset(NMS_labelled_points_l, NMS_labelled_points_l$DecidPer >=20 &NMS_labelled_points_l$DecidPer <60 )
CON_mix<-subset(NMS_labelled_points_l, NMS_labelled_points_l$PinePer >=15 &NMS_labelled_points_l$SprucePer >=15 &NMS_labelled_points_l$FirPer >=15 )

plot(x=NMS_labelled_points_l$X1, y=NMS_labelled_points_l$X2, xlim = c(-.7,.7), ylim = c(-.7,.7), xaxt = "n", yaxt = "n", xlab = "", ylab = "", pch = 16)
points(x=PL_dominant$X1, y=PL_dominant$X2,pch = 10, col = "red")
points(x=BL_dominant$X1, y=BL_dominant$X2,pch = 10,  col = "orange")
points(x=SX_dominant$X1, y=SX_dominant$X2,pch = 2,  col = "green")
points(x=DECID_mix$X1, y=DECID_mix$X2,pch = 3, col = "blue", cex = 1.5)
points(x= CON_mix$X1, y = CON_mix$X2, pch = 3, col = "purple", cex = 1.5)
legend(x="bottomleft", legend = c( "PL_dominant", "BL_dominant", 'SX_dominant','DECID_mix' ,'CON_mix'), pch = c(16,10,2,3, 3), col = c("red", "orange", "green", "blue", "purple"))

#what do we have to work with in age class 2 and 5 for the above categories
#lots of options for age class 2 pine dom, can exclude plots with > 20% deciduous
subset(PL_dominant, PL_dominant$tot_stand_age > 20 & PL_dominant$tot_stand_age <40)
# only stands that are decently stocked and not "T" type... lots of choice
subset(PL_dominant, PL_dominant$tot_stand_age > 20 & PL_dominant$tot_stand_age <40 & PL_dominant$stemsha_LIV>700& PL_dominant$sampletype != "T")
PL_Age2<-subset(PL_dominant, PL_dominant$tot_stand_age > 20 & PL_dominant$tot_stand_age <40& PL_dominant$stemsha_LIV>700 & PL_dominant$DecidPer <15 & PL_dominant$sampletype != "T")
PL_Age2<-subset(PL_Age2, PL_Age2$samp_id %in% sample(PL_Age2$samp_id, 8))
#lots of options for age class 4 pine dom
subset(PL_dominant, PL_dominant$tot_stand_age >= 60 & PL_dominant$tot_stand_age <80)
# lots of age class 5 pine dominant stands available that are decently stocked, combo of VRI and G types
subset(PL_dominant, PL_dominant$tot_stand_age >= 80 & PL_dominant$tot_stand_age <=100& PL_dominant$stemsha_LIV>700)
PL_Age5<-subset(PL_dominant, PL_dominant$tot_stand_age >= 80 & PL_dominant$tot_stand_age <=100& PL_dominant$stemsha_LIV>700& PL_dominant$lead_si1 >13)
PL_Age5<-subset(PL_Age5, PL_Age5$samp_id %in% sample(PL_Age5$samp_id, 8))

# 10 spruce dom age class 2 available, only 10 if we exclude "T" sample and 9 if excluding >15% deciduous
subset(SX_dominant, SX_dominant$tot_stand_age > 20 & SX_dominant$tot_stand_age <40 & SX_dominant$stemsha_LIV>700)
SX_Age2<-subset(SX_dominant, SX_dominant$tot_stand_age > 20 & SX_dominant$tot_stand_age <40 & SX_dominant$stemsha_LIV>700& SX_dominant$sampletype != "T"& SX_dominant$DecidPer <15)
# 7 spruce dom age class 4 available, only 6 if we exclude "T" sample
subset(SX_dominant, SX_dominant$tot_stand_age >= 60 & SX_dominant$tot_stand_age <80 & SX_dominant$stemsha_LIV>700)
# 4 spruce dom age class 5 available, 8 if we go from 80-120
subset(SX_dominant, SX_dominant$tot_stand_age >= 80 & SX_dominant$tot_stand_age <=100& SX_dominant$stemsha_LIV>700)
SX_Age5_6<-subset(SX_dominant, SX_dominant$tot_stand_age >= 80 & SX_dominant$tot_stand_age <=120& SX_dominant$stemsha_LIV>700)

# 1 fir dom age class 2 available
subset(BL_dominant, BL_dominant$tot_stand_age > 20 & BL_dominant$tot_stand_age <40 & BL_dominant$stemsha_LIV>700)
# 5 fir dom age class 4 available
subset(BL_dominant, BL_dominant$tot_stand_age >= 60 & BL_dominant$tot_stand_age <80 & BL_dominant$stemsha_LIV>700)
# 5 fir dom age class 5 available
subset(BL_dominant, BL_dominant$tot_stand_age >= 80 & BL_dominant$tot_stand_age <=100& BL_dominant$stemsha_LIV>700)

# 3 CON mix age class 2 available
subset(CON_mix, CON_mix$tot_stand_age > 20 & CON_mix$tot_stand_age <40 & CON_mix$stemsha_LIV>700)
# 3 CON mix age class 4 available
subset(CON_mix, CON_mix$tot_stand_age >= 60 & CON_mix$tot_stand_age <80 & CON_mix$stemsha_LIV>700)
# 8 CON mix  age class 5 available
subset(CON_mix, CON_mix$tot_stand_age >= 80 & CON_mix$tot_stand_age <=100& CON_mix$stemsha_LIV>700)
CON_Age5<-subset(CON_mix, CON_mix$tot_stand_age >= 80 & CON_mix$tot_stand_age <=100& CON_mix$stemsha_LIV>700)

# 9 DECID mix age class 2 available, 8 if we exclude "T" types
subset(DECID_mix, DECID_mix$tot_stand_age > 20 & DECID_mix$tot_stand_age <40 & DECID_mix$stemsha_LIV>700)
DECID_Age2<-subset(DECID_mix, DECID_mix$tot_stand_age > 20 & DECID_mix$tot_stand_age <40 & DECID_mix$stemsha_LIV>700 & DECID_mix$sampletype != "T")
# 5 DECID mix age class 4 available
subset(DECID_mix, DECID_mix$tot_stand_age >= 60 & DECID_mix$tot_stand_age <80 & DECID_mix$stemsha_LIV>700)
# 2 DECID mix  age class 5 available, 3 stands if we go from 80-140
subset(DECID_mix, DECID_mix$tot_stand_age >= 80 & DECID_mix$tot_stand_age <=100& DECID_mix$stemsha_LIV>700)


plot(x=NMS_labelled_points_l$X1, y=NMS_labelled_points_l$X2, xlim = c(-.7,.7), ylim = c(-.7,.7), xaxt = "n", yaxt = "n", xlab = "", ylab = "", pch = 16, col = "transparent")
points(x=PL_Age2$X1, y=PL_Age2$X2,pch = 10, col = "red")
points(x=PL_Age5$X1, y=PL_Age5$X2,pch = 10, col = "red", cex = 1.5)
points(x=SX_Age2$X1, y=SX_Age2$X2,pch = 2,  col = "green")
points(x=SX_Age5_6$X1, y=SX_Age5_6$X2,pch = 2,  col = "green", cex = 1.5)
points(x=DECID_Age2$X1, y=DECID_Age2$X2,pch = 3, col = "blue")
points(x= CON_Age5$X1, y = CON_Age5$X2, pch = 3, col = "purple", cex = 1.5)
legend(x="bottomleft", legend = c( "PL_dominant", 'SX_dominant','DECID_mix' ,'CON_mix'), pch = c(16,2,3, 3), col = c("red", "green", "blue", "purple"))
#######################################################################################
#########Second ordination with descriptive stand variables included in multivariate matrix
#######################################################################################
SBSmc_tree_output_m<-merge(SBSmc_tree_output, SBSmc_sample[c("samp_id", "tot_stand_age", "dbhq_LIV", 
                                                                "baha_LIV", "stemsha_LIV", "wsvha_LIV", 
                                                                "lead_si1")], by = "samp_id")

#drop samples without values for new variables
SBSmc_tree_output_m<- SBSmc_tree_output_m[is.na(SBSmc_tree_output_m$tot_stand_age) == FALSE,]
min(SBSmc_tree_output_m$dbhq_LIV)
SBSmc_tree_output_m[is.na(SBSmc_tree_output_m$dbhq_LIV),]
SBSmc_tree_output_m<- SBSmc_tree_output_m[is.na(SBSmc_tree_output_m$dbhq_LIV) == FALSE,]
SBSmc_tree_output_m[is.na(SBSmc_tree_output_m$lead_si1),]
SBSmc_tree_output_m<- SBSmc_tree_output_m[is.na(SBSmc_tree_output_m$lead_si1) == FALSE,]
SBSmc_tree_output_m<- subset(SBSmc_tree_output_m, SBSmc_tree_output_m$tot_stand_age >0)

Codes2<-SBSmc_tree_output_m[,1] #save the info for each tree in a codes file

#get bray-curtis distances, run NMDS and output minimum distance configuration with stress and R2 info
dist2<-distance(SBSmc_tree_output_m[,2:length(SBSmc_tree_output_m)], "bray-curtis")
distNMDS2<-nmds(dist2, mindim = 2, maxdim = 2, maxit = 100000)
dist.scors2<-nmds.min(distNMDS2)

#create dataframe with nmds results and attach to codes
NMS_labelled_points2<- data.frame(dist.scors2, Codes2 )
colnames(NMS_labelled_points2) <-c("X1",   "X2", "samp_id")

setwd("C:/Users/elilles/OneDrive - Government of BC/SORTIE_climate_scenarios")
#write.csv(NMS_labelled_points2, file = "NMS_results_PSP_SBSmc_ALL_2.csv")


NMS_labelled_points2_l<-merge(NMS_labelled_points2, SBSmc_sample[c("samp_id", "sampletype", "bgc_ss_grd", "tsa", "meas_yr",  "spc_live1", "tot_stand_age", "dbhq_LIV", "baha_LIV", "stemsha_LIV", "wsvha_LIV", "lead_si1")], by = "samp_id")
NMS_labelled_points2_l<-merge(NMS_labelled_points2_l, SBSmc_tree_output_m[,1:15], by = "samp_id")

young<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$tot_stand_age<= 40)
young_mat<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$tot_stand_age> 40 & NMS_labelled_points2_l$tot_stand_age<= 80)
mature<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$tot_stand_age> 80 &NMS_labelled_points2_l$tot_stand_age<= 140)
old<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$tot_stand_age> 140 & NMS_labelled_points2_l$tot_stand_age<= 230)
v.old<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$tot_stand_age> 230)

plot(x=NMS_labelled_points2_l$X1, y=NMS_labelled_points2_l$X2, xlim = c(-.7,.7), ylim = c(-.7,.7), xaxt = "n", yaxt = "n", xlab = "", ylab = "", pch = 16)
points(x=young$X1, y=young$X2,pch = 10, col = "red")
points(x=young_mat$X1, y=young_mat$X2,pch = 10,  col = "orange")
points(x=mature$X1, y=mature$X2,pch = 2,  col = "green")
points(x=old$X1, y=old$X2,pch = 3, col = "blue")
points(x=v.old$X1, y=v.old$X2,pch = 5, col = "purple")
legend(x="bottomright", legend = c( "<40", "40-80", '80-140','140-230','>230'), pch = c(10,10,2,3,5), col = c("red", "orange", "green", "blue", "purple"))

#add vectors
vectors2=vf(dist.scors2, SBSmc_tree_output_m[,2:length(SBSmc_tree_output_m)], nperm=10)
#vectors<-subset(vectors, vectors[,4] <= 0.1) trying to only plot vectors with p-values of 0.1 or lower
plot(vectors2, 
     #len=0.1,
     col="red")
####

#by datasource
CMI<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$sampletype == "CMI" )
VRI<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$sampletype == "VRI" )
YSM<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$sampletype == "YSM" )
G<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$sampletype == "G" )
Tplots<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$sampletype == "T" )

plot(x=NMS_labelled_points2_l$X1, y=NMS_labelled_points2_l$X2, xlim = c(-.7,.7), ylim = c(-.7,.7), xaxt = "n", yaxt = "n", xlab = "", ylab = "", pch = 16)
points(x=CMI$X1, y=CMI$X2,pch = 10, col = "red")
points(x=VRI$X1, y=VRI$X2,pch = 10,  col = "orange")
points(x=YSM$X1, y=YSM$X2,pch = 2,  col = "green")
points(x=G$X1, y=G$X2,pch = 3, col = "blue")
points(x=Tplots$X1, y=Tplots$X2,pch = 5, col = "purple")
legend(x="bottomright", legend = c( "CMI", "VRI", 'YSM','G','Tplots'), pch = c(16,10,2,3, 5), col = c("red", "orange", "green", "blue", "purple"))



#by leading species
PL_leading<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$spc_live1 == "PL" )
BL_leading<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$spc_live1 == "BL" )
SX_leading<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$spc_live1 == "SX" |NMS_labelled_points2_l$spc_live1 == "SE")
DECID_leading<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$spc_live1 %in% deciduous )
horsetail_SS<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$bgc_ss_grd == 10)

plot(x=NMS_labelled_points2_l$X1, y=NMS_labelled_points2_l$X2, xlim = c(-.7,.7), ylim = c(-.7,.7), xaxt = "n", yaxt = "n", xlab = "", ylab = "", pch = 16)
points(x=PL_leading$X1, y=PL_leading$X2,pch = 10, col = "red")
points(x=BL_leading$X1, y=BL_leading$X2,pch = 10,  col = "orange")
points(x=SX_leading$X1, y=SX_leading$X2,pch = 2,  col = "green")
points(x=DECID_leading$X1, y=DECID_leading$X2,pch = 3, col = "blue")
points(x= horsetail_SS$X1, y = horsetail_SS$X2, pch = 3, col = "purple")
legend(x="bottomright", legend = c( "PL_leading", "BL_leading", 'SX_leading','DECID_leading' ,'horsetail'), pch = c(16,10,2,3, 3), col = c("red", "orange", "green", "blue", "purple"))

#pull out plots that have small tree tallies to go with them
#just G plots that are present in regen database


#NMS_labelled_points2_l.G<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$sampletype == "G")
NMS_labelled_points2_l.G<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$samp_id %in% PSPregen$samp_id)

PL_leading<-subset(NMS_labelled_points2_l.G, NMS_labelled_points2_l.G$spc_live1 == "PL" )
BL_leading<-subset(NMS_labelled_points2_l.G, NMS_labelled_points2_l.G$spc_live1 == "BL" )
SX_leading<-subset(NMS_labelled_points2_l.G, NMS_labelled_points2_l.G$spc_live1 == "SX" |NMS_labelled_points2_l.G$spc_live1 == "SE")
DECID_leading<-subset(NMS_labelled_points2_l.G, NMS_labelled_points2_l.G$spc_live1 %in% deciduous )

plot(x=NMS_labelled_points2_l.G$X1, y=NMS_labelled_points2_l.G$X2, xlim = c(-.7,.7), ylim = c(-.7,.7), xaxt = "n", yaxt = "n", xlab = "", ylab = "", pch = 16)
points(x=PL_leading$X1, y=PL_leading$X2,pch = 10, col = "red")
points(x=BL_leading$X1, y=BL_leading$X2,pch = 10,  col = "orange")
points(x=SX_leading$X1, y=SX_leading$X2,pch = 2,  col = "green")
points(x=DECID_leading$X1, y=DECID_leading$X2,pch = 3, col = "blue")
legend(x="bottomleft", legend = c( "PL_leading", "BL_leading", 'SX_leading','DECID_leading'), pch = c(16,10,2,3), col = c("red", "orange", "green", "blue"))

#################################
##### dominant species
################################
NMS_labelled_points2_l$PinePer <-100*rowSums(NMS_labelled_points2_l[c("PLInit.Dens.10", "PLInit.Dens.20", "PLInit.Dens.30","PLInit.Dens.80")])/  NMS_labelled_points2_l$stemsha_LIV
NMS_labelled_points2_l$SprucePer <-100*rowSums(NMS_labelled_points2_l[c("SXInit.Dens.10", "SXInit.Dens.20", "SXInit.Dens.30","SXInit.Dens.80")])/  NMS_labelled_points2_l$stemsha_LIV
NMS_labelled_points2_l$FirPer <-100*rowSums(NMS_labelled_points2_l[c("BLInit.Dens.10", "BLInit.Dens.20", "BLInit.Dens.30","BLInit.Dens.80")])/  NMS_labelled_points2_l$stemsha_LIV
NMS_labelled_points2_l$DecidPer <-100*rowSums(NMS_labelled_points2_l[c("DECID_sm",  "DECID_lg" )])/  NMS_labelled_points2_l$stemsha_LIV

PL_dominant<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$PinePer >= 70 & NMS_labelled_points2_l$PinePer <=100)
SX_dominant<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$SprucePer >= 70 & NMS_labelled_points2_l$SprucePer <=100 )
BL_dominant<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$FirPer >= 70 & NMS_labelled_points2_l$FirPer <=100)
DECID_mix<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$DecidPer >=20 &NMS_labelled_points2_l$DecidPer <60 )
CON_mix<-subset(NMS_labelled_points2_l, NMS_labelled_points2_l$PinePer >=15 &NMS_labelled_points2_l$SprucePer >=15 &NMS_labelled_points2_l$FirPer >=15 )

plot(x=NMS_labelled_points2_l$X1, y=NMS_labelled_points2_l$X2, xlim = c(-.7,.7), ylim = c(-.7,.7), xaxt = "n", yaxt = "n", xlab = "", ylab = "", pch = 16)
points(x=PL_dominant$X1, y=PL_dominant$X2,pch = 10, col = "red")
points(x=BL_dominant$X1, y=BL_dominant$X2,pch = 10,  col = "orange")
points(x=SX_dominant$X1, y=SX_dominant$X2,pch = 2,  col = "green")
points(x=DECID_mix$X1, y=DECID_mix$X2,pch = 3, col = "blue", cex = 1.5)
points(x= CON_mix$X1, y = CON_mix$X2, pch = 3, col = "purple", cex = 1.5)
legend(x="bottomright", legend = c( "PL_dominant", "BL_dominant", 'SX_dominant','DECID_mix' ,'CON_mix'), pch = c(16,10,2,3, 3), col = c("red", "orange", "green", "blue", "purple"))

#plot chosen stands
PL_2Age2<-subset(PL_dominant, PL_dominant$tot_stand_age > 20 & PL_dominant$tot_stand_age <40& PL_dominant$stemsha_LIV>700 & PL_dominant$DecidPer <15 & PL_dominant$sampletype != "T")
PL_2Age2<-subset(PL_2Age2, PL_2Age2$samp_id %in% sample(PL_2Age2$samp_id, 8))
#add a limitation that site index >13 to get rid of poorish plot 0141_0050_VRI from sample
PL_2Age5<-subset(PL_dominant, PL_dominant$tot_stand_age >= 80 & PL_dominant$tot_stand_age <=100& PL_dominant$stemsha_LIV>700 & PL_dominant$lead_si1 >13)
PL_2Age5<-subset(PL_2Age5, PL_2Age5$samp_id %in% sample(PL_2Age5$samp_id, 8))
SX_2Age2<-subset(SX_dominant, SX_dominant$tot_stand_age > 20 & SX_dominant$tot_stand_age <40 & SX_dominant$stemsha_LIV>700& SX_dominant$sampletype != "T"& SX_dominant$DecidPer <15)
SX_2Age5_6<-subset(SX_dominant, SX_dominant$tot_stand_age >= 80 & SX_dominant$tot_stand_age <=120& SX_dominant$stemsha_LIV>700)
CON_2Age5<-subset(CON_mix, CON_mix$tot_stand_age >= 80 & CON_mix$tot_stand_age <=100& CON_mix$stemsha_LIV>700)
DECID_2Age2<-subset(DECID_mix, DECID_mix$tot_stand_age > 20 & DECID_mix$tot_stand_age <40 & DECID_mix$stemsha_LIV>700 & DECID_mix$sampletype != "T")

plot(x=NMS_labelled_points2_l$X1, y=NMS_labelled_points2_l$X2, xlim = c(-.7,.7), ylim = c(-.7,.7), xaxt = "n", yaxt = "n", xlab = "", ylab = "", pch = 16, col = "transparent")
points(x=PL_2Age2$X1, y=PL_2Age2$X2,pch = 10, col = "red")
points(x=PL_2Age5$X1, y=PL_2Age5$X2,pch = 10, col = "red", cex = 1.5)
points(x=SX_2Age2$X1, y=SX_2Age2$X2,pch = 2,  col = "green")
points(x=SX_2Age5_6$X1, y=SX_2Age5_6$X2,pch = 2,  col = "green", cex = 1.5)
points(x=DECID_2Age2$X1, y=DECID_2Age2$X2,pch = 3, col = "blue")
points(x= CON_2Age5$X1, y = CON_2Age5$X2, pch = 3, col = "purple", cex = 1.5)
legend(x="bottomleft", legend = c( "PL_dominant", 'SX_dominant','DECID_mix' ,'CON_mix'), pch = c(16,2,3, 3), col = c("red", "green", "blue", "purple"))

#final  list of stands
InitiationPlots<-rbind(PL_2Age2, PL_2Age5, SX_2Age2, SX_2Age5_6, DECID_2Age2, CON_2Age5)

#checking for issues... good all the chosen plots have stem sums that check out
InitiationPlots$sumStems<-rowSums(InitiationPlots[,15:28])
plot(InitiationPlots$sumStems, InitiationPlots$stemsha_LIV)
abline(0,1)

#no site series call that would rule plots out
unique(InitiationPlots$bgc_ss_grd)
#check site index
hist(InitiationPlots$lead_si1) #spread is wider that I would like
unique(InitiationPlots$lead_si1) #SI of 12.1 is an issue and 33.7
subset(InitiationPlots, InitiationPlots$lead_si1 <13) #probably OK to keep this spruce dom plot... there are some taller trees and one Pl with Sit index 23, larger dead trees so maybe some site index trees were not open-grown, unlikely to be an '02' site with so much spruce
subset(InitiationPlots, InitiationPlots$lead_si1 >30) #I think we can leave this plot in... SI of sW in the plot is ~24, AT site index is 33.7

#geographic spread of selected plots
InitiationPlots<-merge(InitiationPlots, SBSmc_sample[c("samp_id", "latitude", "longitude", "elev")], by = "samp_id")
hist(InitiationPlots$elev)
InitiationPlots[which(InitiationPlots$elev >1100),]
plot(InitiationPlots$latitude, InitiationPlots$longitude)

SBSmc_tree_output_m_lat_long<-merge(SBSmc_tree_output_m, SBSmc_sample[c("samp_id", "latitude", "longitude", "elev")], by = "samp_id")
plot(SBSmc_tree_output_m_lat_long$latitude, SBSmc_tree_output_m_lat_long$longitude, pch = ".", cex = 2)
points(InitiationPlots$latitude, InitiationPlots$longitude, col = "red", pch =2)

InitiationPlots$StandType<-ifelse(InitiationPlots$samp_id %in% PL_2Age2$samp_id, "PineDom_AgeClass2",
                           ifelse(InitiationPlots$samp_id %in% PL_2Age5$samp_id, "PineDom_AgeClass5",
                           ifelse(InitiationPlots$samp_id %in% SX_2Age2$samp_id, "SpruceDom_AgeClass2",
                           ifelse(InitiationPlots$samp_id %in% SX_2Age5_6$samp_id, "SpruceDom_AgeClass5_6",
                          ifelse(InitiationPlots$samp_id %in% CON_2Age5$samp_id, "ConiferMix_AgeClass5",
                                 "DecidMix_AgeClass2")))))
table(InitiationPlots$StandType, InitiationPlots$sampletype)
write.csv(InitiationPlots, "InitiationPlots_SBSmc.csv")
###some more explanatory plots
NMS_labelled_points2_l$sampletype<-as.factor(NMS_labelled_points2_l$sampletype)
plot(NMS_labelled_points2_l$sampletype, NMS_labelled_points2_l$tot_stand_age, ylab = "Stand2Age")

NMS_labelled_points2_l.G$spc_live1<-as.factor(NMS_labelled_points2_l.G$spc_live1)
plot(NMS_labelled_points2_l.G$spc_live1, NMS_labelled_points2_l.G$tot_stand_age, ylab = "Stand2Age")
plot(NMS_labelled_points2_l.G$spc_live1, NMS_labelled_points2_l.G$lead_si1, ylab = "Site Index")






