#take the psps identified in step 1 and create sortie inits

library(data.table)
library(PSPdataPrep)


in_dir <- file.path("01_sortie_runs","02_parameter_values")
psp_dir <- "D:/Sync/BVRC/FSD program/Data/PSP/PSP data/PSP csv/" #PSP data location
out_dir <- file.path("01_sortie_runs","03_parameter_values")


psp_inits <- fread(file.path(in_dir,"InitiationPlots_SBSmc.csv"))

sample_psp <- import_psps(data_path = psp_dir)

#are the samples all there?
#length(unique(sample_psp[SAMP_ID %in% psp_inits$samp_id]$SAMP_ID))
#get tsas:
selected_psps <- unique(sample_psp[SAMP_ID %in% psp_inits$samp_id]$SAMP_ID)

tsas_sel <- unique(sample_psp[SAMP_ID %in% psp_inits$samp_id]$FileName)

#import trees from tsas and clean
trees_psps <- psp_tree_meas(data_path = psp_dir,
                            tsas = tsas_sel,
                            selected_plots = selected_psps)

#combine plot information and tree data
trees_ss <- merge(trees_psps, psp_inits[,.(samp_id,StandType)], by.x = c("samp_id"),
                  by.y = c("samp_id"), all.x = TRUE, allow.cartesian = TRUE)

#only want last measurement year:
samp_years <- unique(sample_psp[SAMP_ID %in% psp_inits$samp_id & meas_last == "Y",
                                .(SAMP_ID, meas_yr, meas_last)])

trees_ss <- merge(trees_ss, samp_years[,.(SAMP_ID, meas_yr)],
                  by.x = c("samp_id","meas_yr"),
                  by.y = c("SAMP_ID","meas_yr"),
                  all.y = TRUE, allow.cartesian = TRUE)


#3. Clean tree data ------------------------------------------------------------------------------------
unique(trees_ss$species)
trees_ss[,sp_PSP:=ifelse(species=="SW","SX",
                         ifelse(species=="S","SX",
                                ifelse(species=="SE","SX",
                                       ifelse(species=="B","BL",
                                              ifelse(species =="DM","D",
                                                     ifelse(species=="DR","D",
                                                            ifelse(species=="AC","AT",
                                                                   ifelse(species=="E","EP",
                                                                          species))))))))]
unique(trees_ss$sp_PSP)

