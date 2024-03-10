#process outputs

library(foreach)
library(doParallel)
library(data.table)
#---------------------------------------------------------------------------
in_dir <- "./04_sortie_runs/02_outputs/"

files_2_ext <- grep("det",list.files(in_dir, pattern = "SBS",
                                     full.names = FALSE),
                    value = TRUE)

if(!dir.exists(paste0(in_dir,"extracted/"))){
  dir.create(paste0(in_dir,"extracted/"))
}
#if need to extract files:
extractFiles(itype = 1, exname = in_dir, tarnames = files_2_ext)

extFileDir <- paste0(in_dir,"extracted/")
treat_acronym <- "SBS"
treat_parse <- paste0(extFileDir,grep("[[:digit:]].xml$",
                                      grep(treat_acronym,list.files(extFileDir),
                                           value=TRUE),value = TRUE))
plots <- read.csv(file.path("04_sortie_runs","selected_psps.csv"))
plots <- (plots$x)

parallel::detectCores()
numcores <- 38 #76/2
parse_grids <- 0
parse_trees <- 1
# which years to parse
yrs <- seq(0,50) #years vary by plot

t_p <- grep(paste(paste0("_",yrs,".xml"),collapse = "|"),treat_parse, value = TRUE)

#in parallel:--------------

################## make this into a function!!!!!!!!!!!!!!!! ########################

#split treat_parse into treatments for parallel processing
t_p_l <- list()
for(i in 1:length(plots)){
  t_p_l[[i]] <- grep(plots[[i]],t_p, value = TRUE)
}


cl <- parallel::makeCluster(numcores)
doParallel::registerDoParallel(cl)
parallel::clusterEvalQ(cl, c(library(foreach),library(tidyverse),
                             library(data.table),library(rsortie),
                             library(stringr)))
parallel::clusterExport(cl=cl, varlist=c("parse_grids","parse_trees",
                                         "treat_acronym","extFileDir",
                                         "t_p_l"))

g_dt_all <- foreach::foreach(i=1:38)%dopar%{
  #g_dt_all <- foreach::foreach(i=1:1)%dopar%{
  g_dt <- data.table()
  t_dt <- data.table()

  for(ii in 1:length(t_p_l[[i]])){
    # identify which treatment, year and unit is being parsed
    unn <- plots[stringr::str_detect(t_p_l[[i]][ii],plots)]
    yr <- sub('\\.xml$', '',stringr::str_split(t_p_l[[i]][ii],"det_")[[1]][2])
    print(paste("parsing:",unn,"timestep",yr))

    if(parse_grids == 1){
      # parse the output xml grid data
      g <- as.data.table(parseMap(t_p_l[[i]][ii]))

      g[, ':='(timestep = yr, Unit = unn)]
      g_dt <- rbind(g_dt, g, fill=TRUE)
    }

    if(parse_trees == 1){
      # parse the output xml grid data
      t <- as.data.table(parseXML(t_p_l[[i]][ii]))

      t[, ':='(timestep = yr, Unit = unn)]
      t_dt <- rbind(t_dt, t, fill=TRUE)
    }

  }
  if(parse_grids == 1){
    fwrite(g_dt, paste0(extFileDir,treat_acronym,"-",unn,"-grids.csv"), append=FALSE)
  }
  if(parse_trees == 1){
    fwrite(t_dt, paste0(extFileDir,treat_acronym,"-",unn,"-trees.csv"), append=FALSE)
  }
}

parallel::stopCluster(cl)
