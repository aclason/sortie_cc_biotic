library(data.table)


in_dir <- file.path("03_sortie_runs","02_parameter_values")
base_dir <- file.path("03_sortie_runs","01_parameter_files","01_base_files")

out_dir <- file.path("03_sortie_runs","02_parameter_values")


scenarios <- fread(file.path(in_dir,"sbs_scenes_more_stand_types.csv"))
NCI.species <- c("Interior_Spruce","Lodgepole_Pine","Subalpine_Fir","Trembling_Aspen")

#just start with the sbs and climate scenarios --------------------------

#read in the base file:
res <- xml2::read_xml(file.path(base_dir,"SBS.xml"))
xml2::write_xml(res, file.path(base_dir,"temp.xml"))
tmp <- readLines(file.path(base_dir,"temp.xml"), encoding="UTF-8")
xml1 <- gsub("\\\\", "//",tmp)

#values
adult_gr_lines <- rsortie::findFileLine(xml1,
                                             itype = 4,
                                             vargroup = "gr_nciMaxPotentialGrowth",
                                             varname = "gr_nmpgVal species",
                                             varmaster = "NCIMasterGrowth6")

par_val <- c()
for(ii in 1:length(adult_gr_lines)){
  st_start <- stringr::str_locate(xml1[adult_gr_lines[ii]],">")
  st_end <- stringr::str_locate(xml1[adult_gr_lines[ii]],"</")

  par_val[ii] <- as.numeric(substr(xml1[adult_gr_lines[ii]],
                               st_start[1]+1,
                               st_end[1]-1))
}


#change growth - just adult
clim_scenes <- c(0,unique(scenarios$clim_growth))

NCI.species <- c("Interior_Spruce","Lodgepole_Pine","Subalpine_Fir","Trembling_Aspen")

##### Adult NCI Parameter values #####
for(i in 1:length(clim_scenes)){

  # max potential growth - climate influence --------------
  NCI.MaxGrowth <- par_val + (par_val*(clim_scenes[i]/100))

  NCIMasterGrowth6 <- NA
  NCI.adults <- rbind(NCIMasterGrowth6,NCI.MaxGrowth)

  colnames(NCI.adults) <- NCI.species

  write.csv(NCI.adults, paste0(out_dir,"/01_nci_",clim_scenes[i],".csv"), quote=TRUE)
}

