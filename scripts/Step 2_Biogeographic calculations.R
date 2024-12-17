# ------------------------------------------------------------------------------
# Script Name: Step 2_Biogeographic Calculations
# Author: Billy Geary
# Date Created: 17/12/2024
# Last Modified: 17/12/2024
# Purpose: This script uses the DEECA HDMs for all fauna species to undertake a set of biogeographic calculations
#          based on the species range. This includes the number and size of habitat patches and the breadth of pyrome regions
#          that are contained within the species' range. 
# Outputs: A set of csv files with the output of the calculations for each species. Species are ID'd using the VBA Taxon Codes.  
# ------------------------------------------------------------------------------
library(terra)
library(sf)
library(tidyverse)
library(exactextractr)

#####################
#### Pyroregions ####
#####################
pyrome = rast("data_raw/spatial/pyroregion/Australias_pyroregions_raster.tiff")
hdm.files = list.files("/Volumes/BioFutures1/SMP_Data/SMPv4.0/HDMs/Fauna", full=TRUE)

ras = rast(hdm.files[1])
pyrome_proj = project(pyrome, ras, method = "near")
plot(pyrome_proj)

unique_count_fun <- function(x) {
  length(unique(x, na.rm=TRUE))
}

data.out = data.frame()
for(i in hdm.files){
  hdm = rast(i)
  hdm <- ifel(hdm>0, 1, NA)
  dom_pyrome = zonal(pyrome_proj,hdm,fun = modal, na.rm=TRUE)
  breadth_pyrome = zonal(pyrome_proj,hdm,fun = unique_count_fun)
  out = data.frame(filename = names(hdm),
                   dominant_pyrome = dom_pyrome$clusters,
                   pyrome_breadth = breadth_pyrome$clusters)
  data.out = rbind(data.out, out)
  print(i)
}

data.out$TAXON_ID <- str_extract(data.out$filename, "(?<=Spp)\\d+(?=_)")
data.out = data.out %>% select(-filename)

write.csv(data.out, "data_clean/pyrome_species.csv")

###########################
#### EFG Mapping ####
###########################
hdm.files = list.files("/Volumes/BioFutures1/SMP_Data/SMPv4.0/HDMs/Fauna", full=TRUE)
efg =  rast("data_raw/DEECA data/EFG/EFG_NUM_225.tif")

unique_count_fun <- function(x) {
  length(unique(x, na.rm=TRUE))
}

data.out = data.frame()
for(i in hdm.files){
  hdm = rast(i)
  hdm <- ifel(hdm>0, 1, NA)
  dom_efg = zonal(efg,hdm,fun = modal, na.rm=TRUE)
  breadth_efg = zonal(efg,hdm,fun = unique_count_fun)
  out = data.frame(filename = names(hdm),
                   dominant_efg = dom_efg$clusters,
                   efg_breadth = breadth_efg$clusters)
  data.out = rbind(data.out, out)
  print(i)
}

data.out$TAXON_ID <- str_extract(data.out$filename, "(?<=Spp)\\d+(?=_)")
data.out = data.out %>% select(-filename)

write.csv(data.out, "data_clean/efg_species.csv")


###########################
#### Landscape Metrics ####
###########################
library(landscapemetrics)
hdm.files = list.files("/Volumes/BioFutures1/SMP_Data/SMPv4.0/HDMs/Fauna", full=TRUE)

data.out = data.frame()
for(i in hdm.files){
  hdm = rast(i)
  hdm <- ifel(hdm>0, 1, NA)
  # Size of biggest patch
  biggest_patch_size = max(lsm_p_area(hdm)$value, na.rm=TRUE)
  # Number of patches > 10000 Ha
  n_patches = sum(lsm_p_area(hdm)$value > 10000)
  # Nearest neighbour euclidian distance
  # Small values = close average nearest neighbors = less isolated
  # High values = far average nearest neighbors = more isolated
  nned = lsm_l_enn_mn(hdm)$value
  
  out = data.frame(filename = names(hdm),
                        biggest_patch_size = biggest_patch_size,
                        n_patches = n_patches,
                        nned = nned)
  data.out = rbind(data.out, out)
  print(i)
}

write.csv(data.out, "data_clean/patchmetrics_species.csv")
