# ------------------------------------------------------------------------------
# Script Name: Step 4_Compile trait data
# Author: Billy Geary
# Date Created: 17/12/2024
# Last Modified: 17/12/2024
# Purpose: This script compiles trait information for each species in the list generated in Step 1. Traits are compiled from
#          A wide range of sources (see below), including some derived variable sfrom Step 2 and 3. 
# Outputs: A set of csv files with the output of the trait data for each taxon group, as well as a compiled dataset containing data for all groups.  
# ------------------------------------------------------------------------------

#### Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(readxl)
# Species
fauna = read.csv("data_clean/species_list_fauna.csv")

# Trait data
amniote.traits = read.csv("data_raw/Amniote_Database_Aug_2015.csv")
combine.traits = read.csv("data_raw/mammals/combine_trait_data_imputed.csv")
birds.traits = read.csv("data_raw/birds/AVONET1_BirdLife.csv")
frog.traits = read.csv("data_raw/amphibians/AmphiBIO_v1.csv")
# ReptTraits : https://www.nature.com/articles/s41597-024-03079-5#Abs1
reptile.traits =  read_excel("data_raw/reptiles/ReptTraits%20dataset%20v1-1.xlsx")
squambase.traits =  read_excel("data_raw/reptiles/SquamBase1.xlsx")

vicfrog.traits = read_excel("data_raw/amphibians/traits_frogsofvictoria.xlsx")

pyromes = read.csv("data_clean/pyrome_species.csv") %>% select(-X)
fame = read.csv("data_clean/fame_combined_slopes.csv") %>% select(-X, -lm_slope, -ee_lm_slope)

patch = read.csv("data_clean/patchmetrics_species.csv") %>% select(-X)
patch$TAXON_ID <- as.numeric(str_extract(patch$filename, "(?<=Spp)\\d+(?=_)"))
patch = patch %>% select(-filename)

# SMP Traits
load("data_raw/DEECA data/SMP_traits/Mammal_traits.Rds")
load("data_raw/DEECA data/SMP_traits/Bird_traits.Rds")
load("data_raw/DEECA data/SMP_traits/Reptile_traits.Rds")
load("data_raw/DEECA data/SMP_traits/Amphibians_traits.Rds")

# SMP EE Data 
freq.fire = readRDS("data_raw/DEECA data/SMP_EE_planBurn.Rds") # DEECA EE Round 2
ee.data = freq.fire %>% 
  filter(action == "FrequentFRB") %>%
  select(TaxonCode, meanBen) %>%
  group_by(TaxonCode) %>%
  summarise(meanBen_FireFreq = mean(meanBen))

sev.fire = readRDS("data_raw/DEECA data/SMP_EE_futureFire.Rds")
sev.fire = sev.fire %>% distinct(species, TaxonCode, scenarioID, lat, long, refmean, FutureFire)
sev.fire = sev.fire %>% 
  group_by(species, TaxonCode, FutureFire) %>% 
  summarise(meanDN = mean(refmean)) %>%
  pivot_wider(id_cols=c(species, TaxonCode), 
              names_from = FutureFire, values_from = meanDN) %>%
  mutate(FutureImpact = `0` - `3`) %>% drop_na(FutureImpact) %>% ungroup() %>%
  select(TaxonCode, FutureImpact)


################# 
#### Mammals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

mammals <- fauna %>% filter(TaxaGroup=="Mammals")
## Traits to get for each
# Body Mass, Home Range Size, Dispersal
# Habitat preference, Habitat breadth, Stratum, 
# Torpor, Sociality, Burrowing
# Diet Preference, Diet Breadth
# Minimum Number of Young Per Year, Breeding Events Per Year, Minimum Age of Reproduction, Lifespan, Breeding Season
mammal.traits = combine.traits %>% 
  mutate(Binomial = paste(genus, species)) %>%
  select(Binomial, 
         adult_mass_g, home_range_km2, dispersal_km, 
         habitat_breadth_n, foraging_stratum, terrestrial_volant,
         hibernation_torpor, fossoriality, 
         det_diet_breadth_n, dphy_invertebrate, dphy_vertebrate, dphy_plant, 
         freshwater,
         litter_size_n, litters_per_year_n, max_longevity_d) %>% 
  rename("Mass_g" = "adult_mass_g") %>% 
  mutate(n_offspring_year = litter_size_n * litters_per_year_n,
         fossoriality = ifelse(fossoriality > 1, 0, fossoriality))

mammals = left_join(mammals, mammal.traits, by=c("Sci_Name_Alt" = "Binomial"))

mammals <- mammals %>%
  mutate(foraging_stratum = case_when(foraging_stratum == "Ar" ~ "Arboreal",
                                      foraging_stratum == "A" ~ "Aerial",
                                      foraging_stratum == "G" ~ "Terrestrial",
                                      TRUE ~ NA_character_)) %>%
  mutate(Arboreal_Insessorial = case_when(is.na(foraging_stratum) ~ NA_real_,
                                          foraging_stratum == "Arboreal" ~ 1,TRUE ~ 0),
         Aerial = case_when(is.na(foraging_stratum) ~ NA_real_,
                            foraging_stratum == "Aerial" ~ 1,TRUE ~ 0),
         Terrestrial = case_when(is.na(foraging_stratum) ~ NA_real_,
                                 foraging_stratum == "Terrestrial" ~ 1,TRUE ~ 0),
         Aquatic = case_when(is.na(foraging_stratum) ~ NA_real_,
                             freshwater==1 | COMMON_NAME=="Water Rat" ~ 1, TRUE~0)) 

smp.mammals = Mammals_traits %>%
  mutate(diet = case_when(Diet=="C"~"Carnivore",
                          Diet=="H"~"Herbivore",
                          Diet=="I"~"Invertivore",
                          Diet=="INF"~"Infloresence",
                          Diet=="O"~"Omnivore",TRUE~NA),
         diet_carnivore = case_when(is.na(Diet) ~ NA_real_,Diet=="C" ~ 1, TRUE~0),
         diet_herbivore = case_when(is.na(Diet) ~ NA_real_,Diet=="H" ~ 1, TRUE~0),
         diet_invertivore = case_when(is.na(Diet) ~ NA_real_,Diet=="I" ~ 1, TRUE~0),
         diet_infloresence = case_when(is.na(Diet) ~ NA_real_,Diet=="INF" ~ 1, TRUE~0),
         diet_omnivore = case_when(is.na(Diet) ~ NA_real_,Diet=="O" ~ 1, TRUE~0),
         diet_granivore = case_when(is.na(Diet) ~ NA_real_,Diet=="S" ~ 1, TRUE~0)) %>%
  mutate(nesting = nesting.code, 
    nesting = case_when(nesting=="G"~"Ground",
                             nesting=="H"~"Hollows",
                             nesting %in% c("N", "P", "V")~"Branch",
                             nesting=="B"~"Burrow",
                             nesting=="C"~"Cave", TRUE~NA)) %>%
  mutate(nest_burrow = case_when(is.na(nesting) ~ NA_real_,nesting=="Burrow"~1, TRUE~0), 
         nest_cave = case_when(is.na(nesting) ~ NA_real_,nesting=="Cave"~1, TRUE~0), 
         nest_ground = case_when(is.na(nesting) ~ NA_real_,nesting=="Ground"~1, TRUE~0), 
         nest_hollows = case_when(is.na(nesting) ~ NA_real_,nesting=="Hollows"~1, TRUE~0), 
         nest_branch = case_when(is.na(nesting) ~ NA_real_,nesting=="Branch"~1, TRUE~0))
  

mammals = left_join(mammals, smp.mammals, by = c("TAXON_ID"="TaxonCode"))
mammals = left_join(mammals, patch, by = c("TAXON_ID"))
mammals = left_join(mammals, pyromes, by = "TAXON_ID")
mammals = left_join(mammals, ee.data, by = c("TAXON_ID"="TaxonCode"))
mammals = left_join(mammals, sev.fire, by = c("TAXON_ID"="TaxonCode"))
mammals = left_join(mammals, fame, by = c("TAXON_ID"))

# Fix up the bent wing bat southern sub species
mammals[mammals$TAXON_ID=="61343", 44:56] <- mammals[mammals$TAXON_ID=="61342", 44:56]
mammals[mammals$TAXON_ID=="61332", 15] <- 1
vic_mammal_traits = mammals %>%
  mutate(stratum_generalist = 0,
         stratum_aquatic = Aquatic, # Fix this - platypus, water rats 
         stratum_saxicolous = 0, 
         stratum_fossorial = fossoriality, ## Fix this!!! 
         stratum_cryptic = 0,
         Genus = sub("^(\\w+).*", "\\1", Sci_Name)) %>%
        # Taxonomic Information
  select(Taxon_ID = TAXON_ID,
         Scientific_Name = Sci_Name,
         Scientific_Name_Alt = Sci_Name_Alt,
         Genus = Genus,
         Common_Name = COMMON_NAME,
         Taxa_Group = TaxaGroup, 
         Conservation_Listing = Extinction_Risk,
         FFG_Status = FFG_Status,
         # Morphology & Movement
         Mass_g = Mass_g, 
         home_range_km2 = home_range_km2, 
         dispersal_km = dispersal_km,
         
         # Stratum
         stratum = foraging_stratum,
         stratum_aerial = Aerial,
         stratum_arboreal_insessorial = Arboreal_Insessorial,
         stratum_aquatic = stratum_aquatic, 
         stratum_cryptic = stratum_cryptic,
         stratum_fossorial = stratum_fossorial,
         stratum_generalist = stratum_generalist,
         stratum_saxicolous = stratum_saxicolous, 
         stratum_terrestrial = Terrestrial,
         # Nesting
         nesting = nesting, 
         nest_burrow = nest_burrow, 
         nest_cave = nest_cave, 
         nest_ground = nest_ground, 
         nest_hollows = nest_hollows, 
         nest_branch = nest_branch, 
         # Behaviour
         volant = terrestrial_volant, 
         hibernation_torpor = hibernation_torpor, 
         # Diet
         diet = diet,
         diet_breadth_n = det_diet_breadth_n, 
         diet_carnivore = diet_carnivore, 
         diet_herbivore = diet_herbivore,
         diet_invertivore = diet_invertivore, 
         diet_infloresence = diet_infloresence,
         diet_omnivore = diet_omnivore,
         diet_granivore = diet_granivore,
         # Reproduction
         litter_size_n = litter_size_n, 
         litters_per_year_n = litters_per_year_n, 
         max_longevity_d = max_longevity_d, 
         n_offspring_year = n_offspring_year, 
         # Vulnerability to predation
         
         # Vulnerability to Herbivores
         
         # Biogeography
         biggest_patch_size = biggest_patch_size, 
         n_habitat_patches = n_patches, 
         patch_isolation = nned, 
         dominant_pyrome = dominant_pyrome, 
         pyrome_breadth = pyrome_breadth,
         
         # Fire vulnerability
         meanben_firefreq = meanBen_FireFreq, 
         future_fire_impact = FutureImpact, 
         #fame_curve_shape = curve_shape, 
         fame_lm_slope = fame_lm_slope
         )


################# 
#### Birds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
birds <- fauna %>% filter(TaxaGroup=="Birds")
# Birdlife dataset the most complete
birds.traits = birds.traits %>% 
  select(Species1, 
         Mass, 
         Habitat, 
         Primary.Lifestyle,
         Trophic.Niche) %>%
  rename("Mass_g" = "Mass")
birds = left_join(birds, birds.traits, by = c("Sci_Name_Alt" = "Species1"))

amniote.traits.birds = amniote.traits %>%
  mutate(Binomial = paste(genus, species)) %>%
  select(Binomial, 
         litter_or_clutch_size_n,
         maximum_longevity_y, 
         litters_or_clutches_per_y)

birds = left_join(birds, amniote.traits.birds, by = c("Sci_Name_Alt" = "Binomial"))

birds = birds %>% mutate(across(14:16, ~ ifelse(.==-999, NA, .))) %>%
  mutate(Omnivore = case_when(is.na(Trophic.Niche) ~ NA_real_,Trophic.Niche == "Omnivore" ~ 1,TRUE ~ 0),
         Vertivore = case_when(is.na(Trophic.Niche) ~ NA_real_,Trophic.Niche== "Vertivore" ~ 1,TRUE ~ 0),
         Invertivore = case_when(is.na(Trophic.Niche) ~ NA_real_,Trophic.Niche== "Invertivore" ~ 1,TRUE ~ 0),
         Granivore = case_when(is.na(Trophic.Niche) ~ NA_real_,Trophic.Niche== "Granivore" ~ 1,TRUE ~ 0),
         Herbivore_terrestrial = case_when(is.na(Trophic.Niche) ~ NA_real_,Trophic.Niche== "Herbivore terrestrial" ~ 1,TRUE ~ 0),
         Aquatic_predator = case_when(is.na(Trophic.Niche) ~ NA_real_,Trophic.Niche== "Aquatic predator" ~ 1,TRUE ~ 0),
         Nectarivore = case_when(is.na(Trophic.Niche) ~ NA_real_,Trophic.Niche== "Nectarivore" ~ 1,TRUE ~ 0),
         Frugivore = case_when(is.na(Trophic.Niche) ~ NA_real_,Trophic.Niche== "Frugivore" ~ 1,TRUE ~ 0),
         Herbivore_aquatic = case_when(is.na(Trophic.Niche) ~ NA_real_,Trophic.Niche == "Herbivore aquatic" ~ 1,TRUE ~ 0),
         Arboreal_Insessorial = case_when(is.na(Primary.Lifestyle) ~ NA_real_,Primary.Lifestyle== "Insessorial" ~ 1,TRUE ~ 0),
         Terrestrial = case_when(is.na(Primary.Lifestyle) ~ NA_real_,Primary.Lifestyle== "Terrestrial" ~ 1,TRUE ~ 0),
         Generalist = case_when(is.na(Primary.Lifestyle) ~ NA_real_,Primary.Lifestyle== "Generalist" ~ 1,TRUE ~ 0),
         Aerial = case_when(is.na(Primary.Lifestyle) ~ NA_real_,Primary.Lifestyle== "Aerial" ~ 1,TRUE ~ 0),
         Aquatic = case_when(is.na(Primary.Lifestyle) ~ NA_real_,Primary.Lifestyle== "Aquatic" ~ 1,TRUE ~ 0)) 

smp.birds = Birds_traits %>%
  mutate(nesting = NestGuildVicForest) %>%
  mutate(nesting = case_when(nesting=="V" ~ "Branch",
                   nesting=="SH" ~ "Hollows",
                   nesting=="X" & FeedGuild=="W"  ~ "Ground",
                   nesting=='X' & HollowDep==1 ~ "Hollows",
                   nesting=="Ledge" ~ "Ledge", 
                   nesting=="W" ~ "Water", 
                   nesting=="LH" ~ "Hollows", 
                   nesting=="BP" ~ "Brood_Parasite", 
                   nesting=="G" ~ "Ground", 
                   nesting=="B" ~ "Branch", TRUE~NA)) %>%
  mutate(nest_burrow = 0,
         nest_cave = 0,
         nest_ground = case_when(is.na(nesting) ~ NA_real_,nesting=="Ground"~1, TRUE~0), 
         nest_hollows = case_when(is.na(nesting) ~ NA_real_,nesting=="Hollows"~1, TRUE~0), 
         nest_branch = case_when(is.na(nesting) ~ NA_real_,nesting %in% c("Nest", "Branch", "Brood_Parasite","Ledge") ~1, TRUE~0))
 
birds = birds %>% mutate(diet = case_when(Trophic.Niche == "Omnivore" ~ "Omnivore",
                                          Trophic.Niche == "Invertivore" ~ "Invertivore",
                                          Trophic.Niche == "Aquatic predator" ~ "Carnivore",
                                          Trophic.Niche == "Vertivore" ~ "Vertivore", 
                                          Trophic.Niche == "Granviore" ~ "Granivore",
                                          Trophic.Niche == "Nectarivore" ~ "Infloresence", 
                                          Trophic.Niche == "Herbivore aquatic" ~ "Herbivore",
                                          Trophic.Niche == "Herbivore terrestrial" ~ "Herbivore", 
                                          Trophic.Niche == "Frugivore" ~ "Infloresence", TRUE~NA))


birds = left_join(birds, smp.birds, by = c("TAXON_ID"="TaxonCode"))
birds = left_join(birds, patch, by = c("TAXON_ID"))
birds = left_join(birds, pyromes, by = "TAXON_ID")
birds = left_join(birds, ee.data, by = c("TAXON_ID"="TaxonCode"))
birds = left_join(birds, sev.fire, by = c("TAXON_ID"="TaxonCode"))
birds = left_join(birds, fame, by = c("TAXON_ID"))

vic_bird_traits = birds %>%
  mutate(max_longevity_d = maximum_longevity_y*365,
         Genus = sub("^(\\w+).*", "\\1", Sci_Name),
         n_offspring_year = litter_or_clutch_size_n*litters_or_clutches_per_y,
         home_range_km2 = NA,
         dispersal_km = NA,
         volant = 1,
         hibernation_torpor = 0, 
         stratum_aquatic = case_when(HabitatVicForest=="W"~1, TRUE~0), 
         stratum_cryptic = 0, 
         stratum_fossorial = 0,
         stratum_saxicolous = 0,
         nesting = nesting, 
         nest_burrow = nest_burrow, 
         nest_cave = nest_cave, 
         nest_ground = nest_ground, 
         nest_hollows = nest_hollows, 
         nest_branch = nest_branch,
         diet_breadth_n = NA, 
         diet = diet,
         diet_carnivore = max(Vertivore, Aquatic_predator, na.rm=TRUE), 
         diet_herbivore = max(Herbivore_terrestrial, Herbivore_aquatic, na.rm=TRUE),
         diet_invertivore = Invertivore, 
         diet_infloresence = max(Nectarivore, Frugivore, na.rm=TRUE),
         diet_omnivore = Omnivore,
         diet_granivore = Granivore
         
         ) %>%
  # Taxonomic Information
  select(Taxon_ID = TAXON_ID,
         Scientific_Name = Sci_Name,
         Scientific_Name_Alt = Sci_Name_Alt,
         Genus=Genus,
         Common_Name = COMMON_NAME,
         Taxa_Group = TaxaGroup, 
         Conservation_Listing = Extinction_Risk,
         FFG_Status = FFG_Status,
         # Morphology & Movement
         Mass_g = Mass_g, 
         home_range_km2 = home_range_km2, 
         dispersal_km = dispersal_km,
         # Stratum
         stratum = Primary.Lifestyle,
         stratum_aerial = Aerial,
         stratum_arboreal_insessorial = Arboreal_Insessorial,
         stratum_aquatic = stratum_aquatic,
         stratum_cryptic = stratum_cryptic,
         stratum_fossorial = stratum_fossorial,
         stratum_generalist = Generalist,
         stratum_saxicolous = stratum_saxicolous,
         stratum_terrestrial = Terrestrial,
         #Nesting 
         nesting = nesting, 
         nest_burrow = nest_burrow,
         nest_cave = nest_cave,
         nest_ground = nest_ground,
         nest_hollows = nest_hollows,
         nest_branch = nest_branch,
         # Behaviour
         volant = volant, 
         hibernation_torpor = hibernation_torpor, 
         # Diet
         diet = diet,
         diet_breadth_n = diet_breadth_n, 
         diet_carnivore = diet_carnivore,
         diet_herbivore = diet_herbivore,
         diet_invertivore = diet_invertivore,
         diet_infloresence = diet_infloresence,
         diet_omnivore = diet_omnivore,
         diet_granivore = diet_granivore,
         # Reproduction
         litter_size_n = litter_or_clutch_size_n, 
         litters_per_year_n = litters_or_clutches_per_y, 
         max_longevity_d = max_longevity_d, 
         n_offspring_year = n_offspring_year, 
         # Vulnerability to predation
         
         # Vulnerability to Herbivores
         
         # Biogeography
         biggest_patch_size = biggest_patch_size, 
         n_habitat_patches = n_patches, 
         patch_isolation = nned, 
         dominant_pyrome = dominant_pyrome, 
         pyrome_breadth = pyrome_breadth,
         
         # Fire vulnerability
         meanben_firefreq = meanBen_FireFreq, 
         future_fire_impact = FutureImpact, 
         #fame_curve_shape = curve_shape, 
         fame_lm_slope = fame_lm_slope
  )

out = rbind(vic_mammal_traits, vic_bird_traits)
out$nesting = as.factor(out$nesting)
out$stratum = as.factor(out$stratum)
out$diet = as.factor(out$diet)

################# 
#### Reptiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
reptiles <- fauna %>% filter(TaxaGroup=="Reptiles")
reptile.traits = reptile.traits %>% 
  select(Species,
         `Maximum body mass (g)`,
         Microhabitat, `Habitat type`, 
         Diet,
         `Maximum Longevity (years)`, 
         `Mean number of offspring per litter or number of eggs per clutch`, 
         `Number of litters or clutches produced per year`) %>%
  rename("Mass_g" = "Maximum body mass (g)",
         "Habitat" = `Habitat type`,
         maximum_longevity_y = `Maximum Longevity (years)`,
         litter_size_n = `Mean number of offspring per litter or number of eggs per clutch`,
         litters_per_year_n = `Number of litters or clutches produced per year`) %>% 
  mutate(Mass_g = as.numeric(Mass_g))

reptile.traits = reptile.traits %>%
mutate_at(vars("Microhabitat", "Habitat", "Diet", "maximum_longevity_y",
               "litter_size_n", "litters_per_year_n"), ~ na_if(., "NA")) %>%
  mutate_at(vars(maximum_longevity_y, litter_size_n, litters_per_year_n), as.numeric)

reptiles.traits = reptiles %>% 
  select(TAXON_ID, Sci_Name, Sci_Name_Alt) %>% 
  left_join(reptile.traits, by = c("Sci_Name_Alt"="Species"))

squambase.traits = squambase.traits %>% 
  select(`Species name (Binomial)`,
         `maximum mass (derived from allometric equations; in log10, g)`,
         substrate, 
         Diet, `Foraging mode`,
         `Average home range size (minimum; m^2)`,
         `Maximum Longevity (years)`,
         `maximum mean brood size`,
         `Yearly broods (maximum)`)
 
reptiles.traits = left_join(reptiles.traits, squambase.traits, by = c("Sci_Name_Alt" = "Species name (Binomial)"))

reptiles.traits.combined = reptiles.traits %>%
  mutate_at(vars(`Maximum Longevity (years)`, `maximum mean brood size`, `Yearly broods (maximum)`),as.numeric) %>%
  mutate(Microhabitat = coalesce(substrate, Microhabitat),
         Diet = coalesce(Diet.x, Diet.y),
         maximum_longevity_y = coalesce(`Maximum Longevity (years)`, maximum_longevity_y),
         litter_size_n = coalesce(`maximum mean brood size`, litter_size_n),
         litters_per_year_n = coalesce(`Yearly broods (maximum)`, litters_per_year_n)) %>%
  select(TAXON_ID, Sci_Name, Sci_Name_Alt, Mass_g, Microhabitat, Habitat, Diet, maximum_longevity_y, 
         litter_size_n, litters_per_year_n)

reptiles = reptiles %>% 
  left_join(reptiles.traits.combined, by=c("TAXON_ID", "Sci_Name", "Sci_Name_Alt")) %>%
  mutate(Terrestrial = case_when(is.na(Microhabitat) ~ NA_real_,
                                 Microhabitat %in% c("Terrestrial",
                                                     "Arboreal&Terrestrial",
                                                     "Fossorial&Saxicolous&Terrestrial",
                                                     "Fossorial&Terrestrial",
                                                     "Cryptic&Terrestrial",
                                                     "Arboreal&Saxicolous&Terrestrial",
                                                     "Saxicolous&Terrestrial") ~ 1,TRUE ~ 0),
         Arboreal_Insessorial = case_when(is.na(Microhabitat) ~ NA_real_,
                                          Microhabitat %in% c("Arboreal&Saxicolous",
                                                              "Arboreal&Terrestrial",
                                                              "Arboreal&Saxicolous&Terrestrial")~ 1,TRUE ~ 0),
         Saxicolous = case_when(is.na(Microhabitat) ~ NA_real_,
                                Microhabitat %in% c("Arboreal&Saxicolous",
                                                    "Saxicolous&Terrestrial",
                                                    "Fossorial&Saxicolous&Terrestrial",
                                                    "Saxicolous",
                                                    "Arboreal&Saxicolous&Terrestrial")~ 1,TRUE ~ 0),
         Fossorial = case_when(is.na(Microhabitat) ~ NA_real_, 
                               Microhabitat %in% c("Fossorial",
                                                   "Fossorial&Saxicolous&Terrestrial",
                                                   "Fossorial&Terrestrial")~ 1,TRUE ~ 0),
         Aquatic = case_when(is.na(Microhabitat) ~ NA_real_,
                             Microhabitat %in% c("Semi-Aquatic",
                                                 "Aquatic")~ 1,TRUE ~ 0),
         Cryptic = case_when(is.na(Microhabitat) ~ NA_real_,
                             Microhabitat %in% c("Cryptic&Terrestrial",
                                                 "Cryptic")~ 1,TRUE ~ 0),
         Omnivore = case_when(is.na(Diet) ~ NA_real_,Diet == "Omnivorous" ~ 1,TRUE ~ 0),
         Carnivore = case_when(is.na(Diet) ~ NA_real_,Diet== "Carnivorous" ~ 1,TRUE ~ 0),
         Herbivore = case_when(is.na(Diet) ~ NA_real_,Diet== "Herbivorous" ~ 1,TRUE ~ 0)) %>%
  mutate(diet = case_when(Diet == "Omnivorous"~"Omnivore",
                          Diet =="Carnivorous"~"Carnivore",
                          Diet=="Herbivorous"~"Herbivore", TRUE~NA)) %>%
  select(-c(Diet,Microhabitat))

smp.reptiles = Reptiles_traits %>%
  mutate(home_range_km2 = MaleHomerange_aver, 
         dispersal_km = MaleDispersal_distance_aver, 
         AdultDiet = case_when(AdultDiet==""~NA, TRUE~AdultDiet)) 
         
reptiles = left_join(reptiles, smp.reptiles, by = c("TAXON_ID"="TaxonCode"))
reptiles = left_join(reptiles, patch, by = c("TAXON_ID"))
reptiles = left_join(reptiles, pyromes, by = "TAXON_ID")
reptiles = left_join(reptiles, ee.data, by = c("TAXON_ID"="TaxonCode"))
reptiles = left_join(reptiles, sev.fire, by = c("TAXON_ID"="TaxonCode"))
reptiles = left_join(reptiles, fame, by = c("TAXON_ID"))

vic_reptile_traits = reptiles %>%
  mutate(diet = coalesce(AdultDiet,diet)) %>%
  mutate(Omnivore = case_when(is.na(diet) ~ NA_real_,diet %in% c("Omnivore","Carnivore/Omnivore") ~ 1,TRUE ~ 0),
         Invertivore = case_when(is.na(diet) ~ NA_real_,diet == "Invertivore" ~ 1,TRUE ~ 0),
         Carnivore = case_when(is.na(diet) ~ NA_real_,diet %in% c("Carnivore", "Carnivore/Omnivore") ~ 1,TRUE ~ 0),
         Herbivore = case_when(is.na(diet) ~ NA_real_,diet== "Herbivore" ~ 1,TRUE ~ 0)) %>%
  mutate(max_longevity_d = maximum_longevity_y*365,
         Genus = sub("^(\\w+).*", "\\1", Sci_Name),
         n_offspring_year = litter_size_n*litters_per_year_n,
         home_range_km2 = home_range_km2,
         dispersal_km = dispersal_km,
         volant = 0,
         hibernation_torpor = 1, 
         diet_breadth_n = NA, 
         stratum = NA,
         stratum_aerial = 0,
         stratum_generalist = 0,
         nesting = NA, 
         nest_burrow = NA, 
         nest_cave = 0, 
         nest_ground = 0, 
         nest_hollows = 0, 
         nest_branch = 0,
         diet_breadth_n = NA, 
         diet = diet,
         diet_carnivore = Carnivore, 
         diet_herbivore = Herbivore,
         diet_invertivore = Invertivore, 
         diet_infloresence = Herbivore,
         diet_omnivore = Omnivore,
         diet_granivore = Herbivore
  ) %>%
  # Taxonomic Information
  select(Taxon_ID = TAXON_ID,
         Scientific_Name = Sci_Name,
         Scientific_Name_Alt = Sci_Name_Alt,
         Genus=Genus,
         Common_Name = COMMON_NAME,
         Taxa_Group = TaxaGroup, 
         Conservation_Listing = Extinction_Risk,
         FFG_Status = FFG_Status,
         # Morphology & Movement
         Mass_g = Mass_g, 
         home_range_km2 = home_range_km2, 
         dispersal_km = dispersal_km,
         # Stratum
         stratum = stratum,
         stratum_aerial = stratum_aerial,
         stratum_arboreal_insessorial = Arboreal_Insessorial,
         stratum_aquatic = Aquatic,
         stratum_cryptic = Cryptic,
         stratum_fossorial = Fossorial,
         stratum_generalist = stratum_generalist,
         stratum_saxicolous = Saxicolous,
         stratum_terrestrial = Terrestrial,
         #Nesting 
         nesting = nesting, 
         nest_burrow = nest_burrow,
         nest_cave = nest_cave,
         nest_ground = nest_ground,
         nest_hollows = nest_hollows,
         nest_branch = nest_branch,
         # Behaviour
         volant = volant, 
         hibernation_torpor = hibernation_torpor, 
         # Diet
         diet_breadth_n = diet_breadth_n, 
         diet = diet,
         diet_carnivore = diet_carnivore,
         diet_herbivore = diet_herbivore,
         diet_invertivore = diet_invertivore,
         diet_infloresence = diet_infloresence,
         diet_omnivore = diet_omnivore,
         diet_granivore = diet_granivore,
         # Reproduction
         litter_size_n = litter_size_n, 
         litters_per_year_n = litters_per_year_n, 
         max_longevity_d = max_longevity_d, 
         n_offspring_year = n_offspring_year, 
         
         # Vulnerability to predation
         
         # Vulnerability to Herbivores
         
         # Biogeography
         biggest_patch_size = biggest_patch_size, 
         n_habitat_patches = n_patches, 
         patch_isolation = nned, 
         dominant_pyrome = dominant_pyrome, 
         pyrome_breadth = pyrome_breadth,
         
         # Fire vulnerability
         meanben_firefreq = meanBen_FireFreq, 
         future_fire_impact = FutureImpact, 
         #fame_curve_shape = curve_shape, 
         fame_lm_slope = fame_lm_slope
  )

out = rbind(vic_reptile_traits, vic_bird_traits, vic_mammal_traits)


################# 
#### Frogs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
frogs <- fauna %>% filter(TaxaGroup=="Amphibians")
frog.traits = frog.traits %>% 
  select(Species, Body_mass_g, Age_at_maturity_min_y, Longevity_max_y,
         Fos, Ter, Aqu, Arb, Vert, 
         Arthro,
         Leaves, 
         Flowers, 
         Seeds, 
         Fruits,
         Reproductive_output_y) %>% 
  rename("Mass_g" = "Body_mass_g",
         maximum_longevity_y = "Longevity_max_y",
         Fossorial = Fos, 
         Terrestrial = Ter,
         Aquatic = Aqu, 
         Arboreal_Insessorial = Arb)

frogs = left_join(frogs, frog.traits, by = c("Sci_Name_Alt"="Species"))

smp.frogs = Amphibians_traits %>%
  mutate(
    stratum_aquatic_smp = case_when(AquaticAdult>0 ~ 1, TRUE~0),
    nest_burrow_smp = case_when(is.na(Burrow) ~ NA_real_,Burrow > 0 ~ 1, TRUE~0))

vicfrog.traits = vicfrog.traits %>%
  select(Taxon_ID, Scientific_Name, Scientific_Name_Alt, 
         stratum, stratum_aerial, stratum_arboreal_insessorial, stratum_aquatic, stratum_cryptic, stratum_fossorial, stratum_generalist, stratum_saxicolous, stratum_terrestrial,
         nesting, nest_burrow, nest_cave, nest_ground, nest_hollows, nest_branch,
         diet, diet_breadth_n, diet_carnivore, diet_herbivore, diet_invertivore, diet_infloresence, diet_omnivore, diet_granivore,
         litter_size_n, litters_per_year_n)

frogs = left_join(frogs, vicfrog.traits , by = c("TAXON_ID"="Taxon_ID"))

frogs = left_join(frogs, smp.frogs, by = c("TAXON_ID"="TaxonCode"))
frogs = left_join(frogs, patch, by = c("TAXON_ID"))
frogs = left_join(frogs, pyromes, by = "TAXON_ID")
frogs = left_join(frogs, ee.data, by = c("TAXON_ID"="TaxonCode"))
frogs = left_join(frogs, sev.fire, by = c("TAXON_ID"="TaxonCode"))
frogs = left_join(frogs, fame, by = c("TAXON_ID"))

vic_frog_traits = frogs %>%
  mutate(max_longevity_d = maximum_longevity_y*365,
         Genus = sub("^(\\w+).*", "\\1", Sci_Name),
         n_offspring_year = litters_per_year_n*litter_size_n,
         litter_size_n = litter_size_n, 
         litters_per_year_n = litters_per_year_n, 
         home_range_km2 = NA,
         dispersal_km = NA,
         volant = 0,
         hibernation_torpor = NA, #######
         diet_breadth_n = NA, 
         stratum = NA,
         stratum_saxicolous = stratum_saxicolous, 
         stratum_aerial = stratum_aerial,
         stratum_generalist = stratum_generalist,
         stratum_cryptic = stratum_cryptic,
         nesting = NA, 
         nest_burrow = nest_burrow, 
         nest_cave = nest_cave, 
         nest_ground = nest_ground, 
         nest_hollows = 0, 
         nest_branch=0,
         diet_breadth_n = NA, 
         diet = NA,
         diet_carnivore = diet_carnivore, 
         diet_herbivore = diet_herbivore,
         diet_invertivore = diet_invertivore, 
         diet_infloresence = diet_infloresence,
         diet_omnivore = diet_omnivore,
         diet_granivore = diet_granivore
  ) %>%
  # Taxonomic Information
  select(Taxon_ID = TAXON_ID,
         Scientific_Name = Sci_Name,
         Scientific_Name_Alt = Sci_Name_Alt,
         Genus=Genus,
         Common_Name = COMMON_NAME,
         Taxa_Group = TaxaGroup, 
         Conservation_Listing = Extinction_Risk,
         FFG_Status = FFG_Status,
         # Morphology & Movement
         Mass_g = Mass_g, 
         home_range_km2 = home_range_km2, 
         dispersal_km = dispersal_km,
         # Stratum
         stratum = stratum,
         stratum_aerial = stratum_aerial,
         stratum_arboreal_insessorial = stratum_arboreal_insessorial,
         stratum_aquatic = stratum_aquatic_smp, # Aquatic adult
         stratum_cryptic = stratum_cryptic,
         stratum_fossorial = stratum_fossorial,
         stratum_generalist = stratum_generalist,
         stratum_saxicolous = stratum_saxicolous,
         stratum_terrestrial = stratum_terrestrial,
         #Nesting 
         nesting = nesting, 
         nest_burrow = nest_burrow,
         nest_cave = nest_cave,
         nest_ground = nest_ground,
         nest_hollows = nest_hollows,
         nest_branch = nest_branch,
         # Behaviour
         volant = volant, 
         hibernation_torpor = hibernation_torpor, 
         # Diet
         diet_breadth_n = diet_breadth_n, 
         diet = diet,
         diet_carnivore = diet_carnivore,
         diet_herbivore = diet_herbivore,
         diet_invertivore = diet_invertivore,
         diet_infloresence = diet_infloresence,
         diet_omnivore = diet_omnivore,
         diet_granivore = diet_granivore,
         # Reproduction
         litter_size_n = litter_size_n, 
         litters_per_year_n = litters_per_year_n, 
         max_longevity_d = max_longevity_d, 
         n_offspring_year = n_offspring_year, 
         
         # Vulnerability to predation
         
         # Vulnerability to Herbivores
         
         # Biogeography
         biggest_patch_size = biggest_patch_size, 
         n_habitat_patches = n_patches, 
         patch_isolation = nned, 
         dominant_pyrome = dominant_pyrome, 
         pyrome_breadth = pyrome_breadth,
         
         # Fire vulnerability
         meanben_firefreq = meanBen_FireFreq, 
         future_fire_impact = FutureImpact, 
         #fame_curve_shape = curve_shape, 
         fame_lm_slope = fame_lm_slope
  )

#### Compile into one database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vic_fauna_traits = rbind(vic_mammal_traits, vic_bird_traits, vic_reptile_traits, vic_frog_traits)

saveRDS(vic_mammal_traits, "data_clean/vic_mammal_traits.Rds")
saveRDS(vic_bird_traits, "data_clean/vic_birds_traits.Rds")
saveRDS(vic_reptile_traits, "data_clean/vic_reptile_traits.Rds")
saveRDS(vic_frog_traits, "data_clean/vic_frog_traits.Rds")


summary(vic_fauna_traits)

write.csv(vic_fauna_traits, "data_clean/vic_fauna_traits.csv")
saveRDS(vic_fauna_traits, "data_clean/vic_fauna_traits.Rds")

#### ENDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


