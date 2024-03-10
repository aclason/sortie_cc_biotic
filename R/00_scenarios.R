library(data.table)


#sbs
stands <- c("pine", "spruce", "aspen/mixed")
age <- c("young", "mature")
clim_growth <- c(-5, 5, 10, 20)

path <- c("little_doth","lots_doth")
spr_beetle <- c("no_sb","little_sb","lots_sb")
pine_beetle <- c("no_mpb","little_mpb","lots_mpb")
asp_miner <- c("no_mine","little_mine","lots_mine")

base_scene <- CJ(stands, age, clim_growth) #24 stands with climate growth

#if pine dominant? add dothistroma
pine_dom_scene <- CJ("pine",age, clim_growth, path, pine_beetle)
spruce_dom_scene <- CJ("spruce",age, clim_growth, path, spr_beetle)
mixed_dom_scene <- CJ("aspen/mixed",age, clim_growth, path, asp_miner)

setnames(pine_dom_scene, new=c("stand","age","clim_growth","path","insect"))
setnames(spruce_dom_scene, new=c("stand","age","clim_growth","path","insect"))
setnames(mixed_dom_scene, new=c("stand","age","clim_growth","path","insect"))

sbs_scenes <- rbind(pine_dom_scene,spruce_dom_scene,mixed_dom_scene)
fwrite(sbs_scenes, "sbs_scenes_more_stand_types.csv")


#sbs
stands <- c("plantation-pine","plantation-pine-spruce",
            "multi-sp-conifer-dom","multi-sp-decid-dom")
clim_growth <- c(-5, 5, 10, 20)

path <- c("little_doth","lots_doth")
spr_beetle <- c("no_sb","little_sb","lots_sb")
pine_beetle <- c("no_mpb","little_mpb","lots_mpb")
asp_miner <- c("no_mine","little_mine","lots_mine")

base_scene <- CJ(stands, clim_growth) #16 stands with climate growth

pine_dom_scene <- CJ(stands[stands == "multi-sp-conifer-dom"|
                              stands == "plantation-pine"|
                              stands == "plantation-pine-spruce"],
                     clim_growth, path, pine_beetle, spr_beetle)
spruce_dom_scene <- CJ(stands, clim_growth, path, spr_beetle)
mixed_dom_scene <- CJ(stands, clim_growth, path, asp_miner)
