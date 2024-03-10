
library(data.table)
library(PSPdataPrep)
library(ClimateNAr)
library(terra)

psp_dir <- "D:/Sync/BVRC/FSD program/Data/PSP/PSP data/PSP csv/" #PSP data location
out_dir <- file.path("01_sortie_runs","02_parameter_values")


# get PSPs -----------------------------------------------------------------------------------------
#Parameters:
BECzone <- c("SBS","SBPS")
#BECsubzone <- c("mc1","mc2","mc3","dk")
BECsubzone <- NULL
min_time_meas <- 0
DBH_cutoff <- 4
#import all the sample (plot-level) data and select plot ids
sample_psp <- import_psps(data_path = psp_dir)


selected_psps <- select_psps(samples_data = sample_psp,
                            BECzone = BECzone,
                            BECsubzone = BECsubzone,
                            min_remeasure = min_time_meas)

tsas_sel <- unique(sample_psp[SAMP_ID %in% selected_psps]$FileName)

#bring in the trees
trees_psps <- psp_tree_meas(data_path = psp_dir,
                            tsas = tsas_sel,
                            selected_plots = selected_psps)

samp_years <- unique(sample_psp[SAMP_ID %in% selected_psps & meas_last == "Y",
                                .(SAMP_ID, samp_no, meas_yr, meas_last,no_plots, latitude, longitude,
                                  plot_typ, sampletype, tot_stand_age, bgc_zone,bgc_subzone,beclabel)])
samp_years <- samp_years[,.(age = max(tot_stand_age, na.rm = T)),
                         by = .(SAMP_ID, samp_no, meas_yr, meas_last,no_plots,
                                latitude, longitude, plot_typ, sampletype, beclabel)]

#samp_psp_last <- merge(sample_psp, samp_years[,.(SAMP_ID, meas_yr)],
 #                      by.x = c("SAMP_ID","meas_yr"),
  #                     by.y = c("SAMP_ID","meas_yr"),
   #                    all.y = TRUE, allow.cartesian = TRUE)

#variable vs fixed plots
#table(samp_psp_last$sampletype, samp_psp_last$plot_typ)

trees <- merge(trees_psps, samp_years,
                  by.x = c("samp_id","meas_yr"),
                  by.y = c("SAMP_ID","meas_yr"),
                  all.y = TRUE, allow.cartesian = TRUE)
#get rid of na dbhs - just for aging VRI? and just live trees and over DBH cutoff
trees <- trees[!is.na(dbh) & ld %in% c( "L",  "I" )]

trees[, BA_ha := treeCalcs::calc_BA(DBH = dbh)*phf_tree,
      by = seq_len(nrow(trees))]
unique(trees$species)
#let's assume B should be Bl and S == Sx
trees[, species := ifelse(species == "B","BL",
                          ifelse(species == "S", "SX",
                            ifelse(species == "SW", "SX",
                              ifelse(species == "SE", "SX",
                                ifelse(species == "E","EP",
                                 ifelse(species == "DM", "D",
                                  ifelse(species == "DR","D",
                               species)))))))]

stand_sum <- trees[,.(SPH = (sum(phf_tree))/no_plots,
                      BA_ha = (sum(BA_ha))/no_plots),
                          by = .(samp_id,meas_yr,species, plot_typ,no_plots)]

#1. remove plots with > 10 % non standard trees----------
stand_sp <- c("SX", "PL", "AT", "BL")
standard_sp <- stand_sum[species %in% stand_sp,
          .(SPHnorm = sum(SPH), BAnorm = sum(BA_ha)), by = c("samp_id","plot_typ")]
hist(standard_sp[plot_typ =="F"]$BAnorm)
hist(standard_sp[plot_typ =="V"]$BAnorm)

non_standard_sp <- stand_sum[!species %in% stand_sp,
                         .(SPH_not = sum(SPH), BA_not = sum(BA_ha)), by = c("samp_id","plot_typ")]

merg_st_nonst <- merge(standard_sp, non_standard_sp, by = c("samp_id","plot_typ"),
                       all.x = TRUE)
merg_st_nonst[is.na(SPH_not), SPH_not:= 0][is.na(BA_not), BA_not:= 0]

merg_st_nonst[,`:=`(SPH_prop = SPHnorm/(SPHnorm + SPH_not),
                    BA_prop = BAnorm/(BAnorm + BA_not)),
              by = seq_len(nrow(merg_st_nonst))]

#update selected psps by removing samples with > 10% non-standard trees - SPH or BA
selected_psps <- merg_st_nonst[SPH_prop > 0.1]$samp_id
#selected_psps <- merg_st_nonst[BA_prop > 0.1]$samp_id

#2. remove plots with with no age----------
samp_summary <- samp_years[SAMP_ID %in% selected_psps]
selected_psps <- samp_summary[age > 0]$SAMP_ID


#3. remove plots with no or negative site index and select circum mesic
#this is the point where lots get cut - make sure this is the correct behaviour
samp_summary <- samp_years[SAMP_ID %in% selected_psps]
lead_si <- unique(sample_psp[SAMP_ID %in% selected_psps & meas_last == "Y" & !is.na(lead_si1),
                             .(SAMP_ID, lead_si1, spc_live1)])
#show which plots are circum SBSmc2 mesic (SBSmc2 01 site index is 17.9 for
#Spruce and 18.8 for Pine)
lead_si <- lead_si[!is.na(lead_si1) & lead_si1 > 0]
#within 2 of Pine:
circum_mesic_pl <- lead_si[spc_live1 == "PL" & lead_si1 > 16.7 & lead_si1 <= 20.8]$SAMP_ID
circum_mesic_sx <- lead_si[spc_live1 == "SX" & lead_si1 > 15.8 & lead_si1 <= 19.9]$SAMP_ID
selected_psps <- unique(circum_mesic_pl,circum_mesic_sx)

richer <- lead_si[lead_si1 > 20.8]
poorer <- lead_si[lead_si1 < 16.8]

#might have to do this later, to whittle down in large categories, but not exclude from lean cats
samp_summary <- samp_years[SAMP_ID %in% selected_psps]


#4. calculate leading species --------------------------------------
# % leading by stems:
stand_sum[, total_SPH := sum(SPH), by = samp_id]
stand_sum[, perc_SPH := (SPH / total_SPH) * 100, by = .(samp_id, species)]


# % leading by BA:
stand_sum[, total_BA := sum(BA_ha), by = samp_id]
stand_sum[, perc_BA := (BA_ha / total_BA) * 100, by = .(samp_id, species)]

avail_stands <- stand_sum[samp_id %in% samp_summary$SAMP_ID]

avail_stands[species == "PL" & perc_SPH >= 80]
avail_stands[species == "SX" & perc_BA >= 60]

# single species leading
pine_lead <- avail_stands[species == "PL" & perc_BA >= 80]
spruce_lead <- avail_stands[species == "SX" & perc_BA >= 80]

#conifer mix - played with the numbers - the threshold makes a big diff for age class 2
pine_sub <- avail_stands[species == "PL" & perc_BA >= 25]
pine_spruce <- avail_stands[samp_id %in% pine_sub$samp_id & species == "SX" &
                              perc_BA >= 25]
conifer_mix <- avail_stands[samp_id %in% pine_spruce$samp_id]


#deciduous mix
unique(avail_stands$species)
decid <- c("AC", "AT", "EP", "D", "W")
decid_tot <- avail_stands[species %in% decid, .(decid_per_SPH = sum(perc_SPH),
                                                decid_per_BA =sum(perc_BA)), by = "samp_id"]
decid_mix <- decid_tot[decid_per_BA >= 20]
#avail_stands[samp_id %in% decid_mix$samp_id]

samp_summary[, stand_type := ifelse(SAMP_ID %in% pine_lead$samp_id, "pine_lead",
                              ifelse(SAMP_ID %in% spruce_lead$samp_id, "spruce_lead",
                                ifelse(SAMP_ID %in% conifer_mix$samp_id, "con_mix",
                                  ifelse(SAMP_ID %in% decid_mix$samp_id, "decid_mix",
                                                         NA))))]
samp_summary <- samp_summary[!is.na(stand_type)]

#5. add age ------------------------------------------------------------------
samp_summary[, age_bin := ifelse(age<= 20, 1,
                            ifelse(age > 20 & age <= 40, 2,
                             ifelse(age > 40 & age <= 60, 2,
                               ifelse(age > 60 & age <= 80, 4,
                                 ifelse(age > 80 & age <= 100, 5,
                                   ifelse(age > 100 & age <= 140, 6, 7))))))]

table(samp_summary[age_bin == 2| age_bin == 4,.(age_bin,stand_type)])


table(samp_summary[beclabel =="SBSmc2" & stand_type == "pine_lead" & age_bin == 2|
                     beclabel =="SBSmc2" &  stand_type == "pine_lead" & age_bin == 4,
             .(age_bin)])

unique(samp_summary[stand_type == "spruce_lead" & age_bin == 2|
          stand_type == "spruce_lead" & age_bin == 4]$beclabel)

table(samp_summary[stand_type == "spruce_lead" & age_bin == 2|
                     stand_type == "spruce_lead" & age_bin == 4,
                   .(beclabel,age_bin)])

plot_psps <- merge(unique(sample_psp[,.(SAMP_ID, meas_yr, latitude, longitude)]),
                   samp_summary, by = c("SAMP_ID","meas_yr"), all.y = TRUE)
plot_psps <- plot_psps[age_bin == 2|age_bin == 4]

#6. decide on which ecosystems to cut ----------------------------------------
unique(plot_psps$beclabel)
keep_zones <- c("SBSmc2","SBSdk", "SBSdk2","SBSmc3", "SBPSmc", "SBPSdc",
                "SBSmw","SBSmk1")
plot_psp_bec <- plot_psps[beclabel %in% keep_zones]

table(plot_psp_bec[,.(age_bin,stand_type)])


#7. making maps----------------------------------------------------------------

plot_psps_2 <- plot_psps[age_bin == 2]
plot_psps_4 <- plot_psps[age_bin == 4]

#library(sf)
plot_psps_sf <- sf::st_as_sf(plot_psps, coords = c("longitude","latitude"), crs = 4326)
plot_psps_sf <- st_make_valid(plot_psps_sf)
st_is_longlat(plot_psps_sf)
plot_psps_sf_t <- st_transform(plot_psps_sf, crs = 3005)


plot_psps_sf_2 <- sf::st_as_sf(plot_psps_2, coords = c("longitude","latitude"), crs = 4326)
plot_psps_sf_2 <- st_make_valid(plot_psps_sf_2)
st_is_longlat(plot_psps_sf_2)
plot_psps_sf2_t <- st_transform(plot_psps_sf_2, crs = 3005)

plot_psps_sf_4 <- sf::st_as_sf(plot_psps_4, coords = c("longitude","latitude"), crs = 4326)
plot_psps_sf_4 <- st_make_valid(plot_psps_sf_4)
st_is_longlat(plot_psps_sf_4)
plot_psps_sf4_t <- st_transform(plot_psps_sf_4, crs = 3005)


sf::st_write(plot_psps_sf_t, file.path(out_dir, "plot_psps_sf_t.gpkg"), append = FALSE)

sf::st_write(plot_psps_sf2_t, file.path(out_dir, "plot_psps_sf2_t.gpkg"), append = FALSE)
sf::st_write(plot_psps_sf4_t, file.path(out_dir, "plot_psps_sf4_t.gpkg"), append = FALSE)



site_simp <- unique(sample_psp[SAMP_ID %in% selected_psps,.(SAMP_ID,meas_yr,
                                                            latitude,longitude,elev,
                                                            beclabel, tot_stand_age)])
site_simp[SAMP_ID == "47007 R000503"]

stand_types <- merge(stand_types, site_simp,
                     by.x = c("samp_id","meas_yr"), by.y = c("SAMP_ID","meas_yr"),
                     all.x = TRUE)

#hist(stand_types$tot_stand_age)
#stand_types[samp_id == "47007 R000503"]

ggplot(stand_types[samp_id %in% selected_psps[1:10]])+
  geom_boxplot(aes(x = species, y = SPH, fill = species))+
  facet_wrap("samp_id")

dom_trees <- stand_types[, .SD[which.max(SPH)], by = .(samp_id, meas_yr)]

dom_trees <- dom_trees[,.N, by = .(samp_id, meas_yr, species)]
dom_trees[,sum(N), by = "species"]

ggplot(dom_trees)+
  geom_bar(aes(x = species), stat = "count")
#so a lot of pine dominant stands

# Erica's version: major scripting effort to make this re-useable and functional for the project
#bind three TSA datasets together, use sample with a 4 cm dbh compile limit.. available for all plots, only G and t have 2 cm dbh limit
SBSmc_sample <- rbind(TSA_03_sample, TSA_14_sample, TSA_20_sample, TSA_26_sample, TSA_24_sample)
#use the most recent measurement only
SBSmc_sample <- subset(SBSmc_sample, SBSmc_sample$meas_last == "Y")

