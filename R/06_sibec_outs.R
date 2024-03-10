library(data.table)

in_dir <- file.path("02_fits")
out_dir <- "03_parameter_values"


#SORTIE parameter name	R parameter name	Interior spruce	Lodgpole  pine	Subalpine fir	Trembling aspen
plot_seasonal_precipitation = c()
plot_n_dep = c()

cs_nci_ad <- data.table(plot_precip_mm_yr = c(71,	40,	19,	60),
                        PrecipEffAl = c(0.0023,	0.0030,	0.0492,	0.0092),
                        PrecipEffBl = c(2.89,	2.37,	2.12,	3.62),
                        PrecipEffCl = c(4.43,	1.89,	2.94,	1.45),
                        PrecipEffAh = c(0.06,	0.02,	0.07,	0.10),
                        PrecipEffBh = c(7.90,	7.84,	7.78,	0.23),
                        PrecipEffCh = c(50.00,	11.91, 49.56,	0.01),
                        TempEffAl = c(0.05,	0.39,	0.29,	0.00),
                        TempEffBl = c(1612,	872,	917,	930),
                        TempEffCl = c(5.04,	11.05,	17.95,	5.95),
                        TempEffAh = c(0.86,	0.94,	0.91,	0.71),
                        TempEffBh = c(1035,	838,	1161,	1185),
                        TempEffCh = c(179,	258,	484,	292),
                        NitrogenX0 = c(5.4,	35.6,	26.2,	24.2),
                        NitrogenXb = c(4.0,	18.4,	16.3,	15.9),
                        Std.Dev.Norm.Log.Adj = 	c(1.9,	0.7,	1.3,	1.0))
cs_nci_adult <- data.table(t(cs_nci_ad))

colnames(cs_nci_adult) <- c("Interior_Spruce","Lodgepole_Pine","Subalpine_Fir",
                                      "Trembling_Aspen")
cs_nci_adult <- rbind(NA,cs_nci_adult, fill = TRUE)
cs_nci_adult[, x := NULL]
rownames(cs_nci_adult) <- c("NA",
                            "plot_precip_mm_yr",
                            "PrecipEffAl",
                            "PrecipEffBl",
                            "PrecipEffCl",
                            "PrecipEffAh",
                            "PrecipEffBh",
                            "PrecipEffCh",
                            "TempEffAl",
                            "TempEffBl",
                            "TempEffCl",
                            "TempEffAh",
                            "TempEffBh",
                            "TempEffCh",
                            "NitrogenX0",
                            "NitrogenXb",
                            "Std.Dev.Norm.Log.Adj")


fwrite(cs_nci_adult, file.path(out_dir,"cs_nci_adult.csv"), row.names = TRUE)
