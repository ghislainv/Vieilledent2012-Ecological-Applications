#!/usr/bin/R

## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
## web             :https://ghislainv.github.io
## license         :GPLv3
## ==============================================================================

## Libraries
require(readr)
require(dplyr)
require(here)
require(ggplot2)
require(glue)
require(fuzzyjoin)
    
## Load data
df_lab <- read_delim("DataLab-G.csv", delim="\t")
df_moist <- read_delim("DataMoisture-G.csv", delim="\t")
df_gps <- read_csv("DataSites-G.csv")

## Mean residual moisture after drying
df_res_moist <- df_moist |>
    mutate(FType=if_else(Site=="Fort-Dauphin-Dry", "dry", "moist-wet")) |>
    group_by(FType) |>
    summarise(res_moist=mean(Moisture, na.rm=TRUE))
rmoist_d <- df_res_moist$res_moist[df_res_moist$FType=="dry"]
rmoist_m <- df_res_moist$res_moist[df_res_moist$FType=="moist-wet"]

## Correcting dry weight and compute wood density
df_lab <- df_lab |>
    mutate(ExactDryWeight = DryWeight - DryWeight *
               if_else(Site=="Fort-Dauphin-Dry", rmoist_d, rmoist_m) / 100) |>
    mutate(wsg = ExactDryWeight / MoistVol2)

## Add corresponding locality names
df_gps <- df_gps %>%
    mutate(Site=c("Bealanana", "Ivohibe", "Fort-Dauphin-Moist", "Fort-Dauphin-Dry", "Fandriana"))

## CombiningXS data-sets
df_wsg <- df_lab %>%
    left_join(df_gps, by="Site") |>
    group_by(Site) |>
    mutate(treeID=sprintf("%03d",as.numeric(as.factor(Tree)))) |>
    ungroup() |>
    mutate(treeID=paste(treeID, Forest_site, sep="_")) |>
    mutate(TreePart=if_else(TreePart %in% c("b", "tr.b"), "trunk_base", TreePart)) |>
    mutate(TreePart=if_else(TreePart %in% c("h", "tr.h"), "trunk_top", TreePart)) |>
    mutate(TreePart=if_else(TreePart %in% c("m", "tr.m"), "trunk_middle", TreePart)) |>
    mutate(TreePart=if_else(TreePart=="br", "branch", TreePart)) |>
    mutate(Species=if_else(Species=="Grewia aprima", "Grewia apetala", Species)) |>
    mutate(wsg=round(wsg, 3)) |>
    mutate(year=2010) |>
    dplyr::filter(!is.na(wsg)) |>
    rename(latin_binomial=Species) |>
    select(treeID, Family, latin_binomial, TreePart, wsg, Forest_site, Village, Lat, Long)

## ======================
## Cleaning taxonomy
## ======================

## Species list
df_sp <- read_csv("species_list.csv")
df_sp$latin_binomial <- unlist(lapply(strsplit(df_sp$Species, split=" "), function(x) {paste(x[1], x[2], sep=" ")}))
df_sp <- df_sp |>
    select(Species, latin_binomial, Family) |>
    distinct()

## Fuzzy matching
df2 <- df_wsg |>
    stringdist_join(df_sp, by="latin_binomial", mode="left",
                    method="jw", distance_col="dist", max_dist=99) |>
    group_by(latin_binomial.x) |>
    slice_min(order_by=dist, n=1) |>
    ungroup() |>
    select(treeID, Species, latin_binomial.y, Family.y, TreePart, wsg, Forest_site, Village, Lat, Long) |>
    arrange(Forest_site, treeID, TreePart) |>
    rename(tree=treeID, species=Species, latin_binomial=latin_binomial.y, family=Family.y,
           tree_part=TreePart, region=Forest_site, village=Village, lat=Lat, lon=Long)

## Group by tree
df_tree <- df2 |>
    group_by(tree) |>
    mutate(wsg_mean=round(mean(wsg, na.rm=TRUE), 3)) |>
    ungroup() |> select(-wsg, -tree_part) |> distinct() |>
    rename(wsg=wsg_mean)

## Save results
write_csv(df2, file="data_wsg_localities_mada.csv")
write_csv(df_tree, file="data_wsg_localities_mada_tree.csv")

## ======================
## Some statistics
## ======================

nobs <- length(df2$wsg) # 1882
ntree <- length(unique(df2$tree)) # 478
ntaxon <- length(unique(df2$latin_binomial)) # 78

## End
