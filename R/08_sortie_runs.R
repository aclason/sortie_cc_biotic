library(data.table)
library(rsortie)


param_dir <- "./01_sortie_runs/02_parameter_values/"
base_param_dir <- "./01_sortie_runs/01_parameter_files/01_base_files/"
param_file_dir <- "01_sortie_runs/01_parameter_files/"
Output_path <- "D:\\GitHub\\sortie_cc_biotic\\01_sortie_runs\\03_run_outputs\\"

plots <- paste0(read.csv(file.path(param_dir,"selected_psps.csv"))$x,".csv")

#merge behaviour files into 1 and and output directory:
nci <- purrr::map(list.files(param_dir, "_nci", full.names = TRUE), fread)

nci <- purrr::map(nci, function(x)rbind(x, data.table(V1="ShortOutput",
                                   Interior_Spruce= Output_path),
                                   fill=TRUE))
nci <- purrr::map(nci, function(x)rbind(x, data.table(V1="Output",
                                   Interior_Spruce= Output_path),
                                   fill=TRUE))
nci <- purrr::map(nci, function(x)setnames(x, "V1",""))

pf_names <- list.files(param_dir, "_nci", full.names = TRUE)
purrr::map2(nci, pf_names, ~ write.csv(.x, file = .y, quote=TRUE, row.names = FALSE))



lstOfFiles <- data.frame("type"=c(0,
                                  rep(1,length(plots)),
                                  rep(2,length(nci))),
                         "name"=c("SBS.xml",
                                  plots,
                                  list.files(param_dir, "_nci", full.names = FALSE)))


rsortie::makeFiles(lstFiles = lstOfFiles,
                   path_basexmls = base_param_dir,
                   path_newxmls = param_file_dir,
                   path_newvals = param_dir)

# update the length of run ---------------------------------
# run for 1 time step first
files2run <- list.files(param_file_dir,pattern = "SBS", full.names = TRUE)
updateNumYears(files2run[1:35], 50)

batch_files <- files2run[1:35]
runSortiePar(fname = batch_files, numcores = 35, sortie_loc=0)
