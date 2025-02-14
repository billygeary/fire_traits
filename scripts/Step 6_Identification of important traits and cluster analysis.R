# ------------------------------------------------------------------------------
# Script Name: Step 6_Identification of important traits and clusters
# Author: Billy Geary
# Date Created: 17/12/2024
# Last Modified: 17/12/2024
# Purpose: This script uses Boosted Regression Trees to identify the key traits 
#          in each taxonomic group that explain a large amount of the variation
#          in a species' vulnerability to fire, as measured by two seperate response 
#          variables: 1) SMP Benefit of Reducing Fire Frequency and 2) FAME
#          response to time since fire. BRTs are used to estimate the variable importance 
#          for each trait in each taxon group. The variable importance 
#          estimates are then used in cluster analysis to cluster species into groups based 
#          on their fire response trait syndrome (i.e. species that have
#          similar responses to fire).
# Outputs: A series of plots identifying the relative importance of each trait in explaining
#          each species' response to fire. A csv with the trait information from Step 5, but
#          with extra columns specifying the cluster in which each species has been assigned
#          based on the cluster analyses performed. 
# ------------------------------------------------------------------------------

#################
#### Setup ###### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
library(tidyverse)
library(gbm)
library(brms)
library(ggplot2)
library(pdp)
library(vegan)
library(mclust)
library(cluster)
library(cowplot)

source("scripts/brt.functions.R")
source("scripts/trait_analysis_functions.R")

# Read in the imputed data
vic_fauna_traits_imputed = readRDS("data_clean/trait_outputs/vic_fauna_traits_imputed.Rds")

vic_fauna_traits_imputed = vic_fauna_traits_imputed  %>% 
  mutate(across(c("nesting", "nest_burrow","nest_cave","nest_ground","nest_hollows","nest_branch",
                  "stratum","stratum_aerial","stratum_arboreal_insessorial", "stratum_aquatic",
                  "stratum_cryptic","stratum_fossorial", "stratum_terrestrial", "stratum_saxicolous",
                  "diet","diet_breadth_n","diet_carnivore","diet_herbivore",
                  "diet_invertivore","diet_infloresence","diet_omnivore","diet_granivore",
                  "dominant_pyrome", "volant", "hibernation_torpor"), as.factor))

# What is our data coverage for each variable
summary_data = vic_fauna_traits_imputed %>% 
  group_by(Taxa_Group) %>%
  summarise(recent_fire = sum(!is.na(fame_lm_slope)),
            frequent_fire = sum(!is.na(meanben_firefreq)),
            severe_fire = sum(!is.na(future_fire_impact)),
            total_spp = n())
summary_data
#################
#### Mammals #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
#### Find Important Traits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Part 1 - Select the traits for inclusion
mammal_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Mammals")

# Define the trait groups
trait_groups = list(Movement_Physiology = c("Mass_g", "home_range_km2", "dispersal_km",  "max_longevity_d"),
                    Reproduction = c("litter_size_n", "litters_per_year_n","n_offspring_year"),
                    Broad_Habitat = c("dominant_pyrome"),
                    Microhabitat = c("stratum","stratum_aerial","stratum_arboreal_insessorial", "stratum_aquatic","stratum_cryptic","stratum_fossorial", "stratum_terrestrial", "stratum_saxicolous"),
                    Nesting = c("nesting", "nest_burrow","nest_cave","nest_ground","nest_hollows","nest_branch"),
                    Behaviour = c("hibernation_torpor", "volant"),
                    Diet = c("diet","diet_breadth_n","diet_carnivore","diet_herbivore",
                             "diet_invertivore","diet_infloresence","diet_omnivore","diet_granivore"),
                    Biogeography = c("biggest_patch_size","n_habitat_patches","patch_isolation", "pyrome_breadth"))

# Convert the list to a data frame using stack
trait_df <- stack(trait_groups)

# Rename columns for clarity
colnames(trait_df) <- c("var", "Group")

# Traits of interest
all_traits <- c("Mass_g", "home_range_km2", "dispersal_km", 
                "stratum_aerial", "stratum_arboreal_insessorial", "stratum_aquatic", "stratum_terrestrial", 
                "nest_burrow", "nest_cave", "nest_ground", "nest_hollows", "nest_branch", 
                "diet_carnivore", "diet_herbivore", "diet_invertivore", "diet_infloresence", "diet_omnivore", 
                "volant", "hibernation_torpor", 
                "litter_size_n", "litters_per_year_n", "max_longevity_d", "n_offspring_year",
                "biggest_patch_size", "n_habitat_patches", "dominant_pyrome", "pyrome_breadth")

all_traits_clean_names = data.frame(var = all_traits, 
                                    traits_label = c("Mass (g)", "Home Range (km2)", "Dispersal (km)", 
                                                     "Stratum - Aerial", "Stratum - Arboreal", "Stratum - Aquatic", "Stratum - Terrestrial", 
                                                     "Nesting - Burrow","Nesting - Cave","Nesting - Ground","Nesting - Hollows","Nesting - Branches",
                                                     "Diet - Carnivore", "Diet - Herbivore", "Diet - Invertivore", "Diet - Infloresence", "Diet - Omnivore", 
                                                     "Volant", "Hibernation/Torpor","Litter Size", "Litters/Yr", "Longevity", "Offspring/Yr", 
                                                     "Size of Largest Habitat Patch", "Number of Habitat Patches", "Dominant Pyrome", "Pyrome Breadth"))
trait_df = left_join(trait_df, all_traits_clean_names)

# Create the formula for each model
traits_formula_smp <- as.formula(paste("meanben_firefreq", "~", paste(all_traits, collapse = " + ")))
traits_formula_fame <- as.formula(paste("fame_lm_slope", "~", paste(all_traits, collapse = " + ")))

#### Too Frequent Fire (SMP)
smp_modelling = mammal_traits %>% select(meanben_firefreq, all_of(all_traits)) %>%
  drop_na(meanben_firefreq)
mammal_brt_smp <- gbm.step(data = smp_modelling, 
                           gbm.x = 2:28,
                           gbm.y = 1,
                           family = "gaussian", 
                           tree.complexity = 2, 
                           learning.rate = 0.001,
                           bag.fraction = 0.5)

# Look at the variable importance
mammal.importance = summary(mammal_brt_smp)  # Variable importance
all_traits_clean_names$traits_label = as.factor(all_traits_clean_names$traits_label)

smp.mammal.importance = mammal.importance %>%
  left_join(trait_df, by = "var") %>%
  mutate(traits_label = fct_reorder(traits_label, rel.inf, .desc = FALSE)) %>%  # Reorder `var` based on `rel.inf`
  mutate(fire.var = "Frequent Fire")

# Partial Dependence Plots
factor_vars <- names(smp_modelling)[sapply(smp_modelling, is.factor)==TRUE]
pdp_data <- lapply(mammal_brt_smp$var.names, get_pdp_data, model =mammal_brt_smp)
mammal_pdps <- lapply(pdp_data, create_pdp_plots, factor_vars)
mammal_smp_pdp <- wrap_plots(mammal_pdps, ncol = 4)  # 5 rows, 3 columns
ggsave(plot = mammal_smp_pdp, filename = "plots/mammal.smp.pdp.jpg", width = 8, height = 5, scale=2)

#### Recent Fire (FAME Slopes)
# Fit the model as a Boosted Regression Tree
fame_modelling = mammal_traits %>% select(fame_lm_slope, all_of(all_traits)) %>%
  drop_na(fame_lm_slope) %>% select(-c(stratum_aquatic, nest_cave)) %>%
  mutate(fame_lm_slope = 1000*fame_lm_slope)

mammal_brt_fame <- gbm.step(data = fame_modelling, 
                           gbm.x = 2:26,
                           gbm.y = 1,
                           family = "gaussian", 
                           tree.complexity = 1, 
                           learning.rate = 0.001,
                           bag.fraction = 0.75)

summary(mammal_brt_fame)

# Look at the variable importance
mammal.importance = summary(mammal_brt_fame)  # Variable importance
all_traits_clean_names$traits_label = as.factor(all_traits_clean_names$traits_label)

fame.mammal.importance = mammal.importance %>%
  left_join(trait_df, by = "var") %>%
  mutate(traits_label = fct_reorder(traits_label, rel.inf, .desc = FALSE)) %>%  # Reorder `var` based on `rel.inf`
  mutate(fire.var = "Recent Fire")

# Partial Dependence Plots
factor_vars <- names(fame_modelling)[sapply(fame_modelling, is.factor)==TRUE]
pdp_data <- lapply(mammal_brt_fame$var.names, get_pdp_data, model =mammal_brt_fame)
mammal_pdps <- lapply(pdp_data, create_pdp_plots, factor_vars)
mammal_fame_pdp <- wrap_plots(mammal_pdps, ncol = 4)  # 5 rows, 3 columns
mammal_fame_pdp

ggsave(plot = mammal_fame_pdp, filename = "plots/mammal.fame.pdp.jpg", width = 8, height = 5, scale=2)

# Plot the results
mammal.importance.data = rbind(fame.mammal.importance, smp.mammal.importance)
mammal.importance.plot = mammal.importance.data %>% 
  filter(rel.inf > 0) %>%
  mutate(fire.var = factor(fire.var, levels = c("Recent Fire", "Frequent Fire"))) %>%
  ggplot(aes(x = traits_label, y = rel.inf, fill = Group)) + geom_col() + 
  labs(fill = "Trait Group", x="Trait", y = "Relative Influence (%)", title="Mammal Trait Importance") + facet_wrap(~fire.var, nrow=1) +
  scale_fill_viridis_d() + theme_minimal() +
  coord_flip() 
mammal.importance.plot
ggsave(plot = mammal.importance.plot, filename = "plots/combined.mammal.importance.plot.jpg", width = 8, height = 4)

#### Clustering Species by Important Traits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot NMDS with clusters

#### Too Frequent Fire (SMP)
smp_traits_scaled = mammal_traits %>% 
  select(all_of(all_traits)) %>% 
  select(all_of(smp.mammal.importance$var)) # Order to match with vector of importance
smp.mammal.importance = smp.mammal.importance$rel.inf/100

# Compute Gower distance matrix
smp_gower_dist <- daisy(smp_traits_scaled, metric = "gower", weights = smp.mammal.importance)
smp_gmm_model <- Mclust(smp_gower_dist, G=1:15) # Try between 1-20 clusters to see which is best. 

# Add cluster assignments to the dataset
mammal_traits$smp_cluster <- smp_gmm_model$classification

# Summarize meanben_firefreq by cluster
species_clusters <- mammal_traits %>% select(Common_Name, smp_cluster)

trait_summary <- mammal_traits %>%
  group_by(smp_cluster) %>%
  summarise(mean_firefreq = mean(meanben_firefreq, na.rm = TRUE),
            sd_firefreq = sd(meanben_firefreq, na.rm = TRUE))

ggplot(mammal_traits) + geom_boxplot(aes(x=as.factor(smp_cluster), y = meanben_firefreq))

## NMDS
smp_nmds <- metaMDS(smp_gower_dist,method="gower", k = 2, trymax = 30)

# Extract NMDS scores and perform clustering
smp_nmds_scores <- as.data.frame(scores(smp_nmds))
smp_nmds_scores$Species = mammal_traits$Common_Name
smp_nmds_scores$Cluster = as.factor(mammal_traits$smp_cluster)

#### Recent Fire (FAME Slopes)
fame_traits_scaled = mammal_traits %>% 
  select(all_of(all_traits)) %>% 
  select(all_of(fame.mammal.importance$var)) # Order to match with vector of importance
fame.mammal.importance = fame.mammal.importance$rel.inf/100

# Compute Gower distance matrix
fame_gower_dist <- daisy(fame_traits_scaled, metric = "gower", weights = fame.mammal.importance)
fame_gmm_model <- Mclust(fame_gower_dist, G=1:15) # Try between 1-20 clusters to see which is best. 

# Add cluster assignments to the dataset
mammal_traits$fame_cluster <- fame_gmm_model$classification

# Summarize meanben_firefreq by cluster
species_clusters <- mammal_traits %>% select(Common_Name, fame_cluster)

ggplot(mammal_traits) + geom_boxplot(aes(x=as.factor(fame_cluster), y =fame_lm_slope))

## NMDS
fame_nmds <- metaMDS(fame_gower_dist,method="gower", k = 2, trymax = 30)

# Extract NMDS scores and perform clustering
fame_nmds_scores <- as.data.frame(scores(fame_nmds))
fame_nmds_scores$Species = mammal_traits$Common_Name
fame_nmds_scores$Cluster = as.factor(mammal_traits$fame_cluster)

# Plot the things
(smp.mammal.pcoa = plot_pcoa(smp_nmds_scores))
(fame.mammal.pcoa = plot_pcoa(fame_nmds_scores))
(mammal.pcoas = plot_grid(fame.mammal.pcoa, smp.mammal.pcoa, labels=c("a)", "b)")))
ggsave(plot = mammal.pcoas, filename = "plots/mammal.pcoa.plot.jpg", width = 8, height = 3.5, scale =1.5)

## Add the traits to the plot
# Compute correlations of traits with NMDS axes
smp_trait_correlations <- envfit(smp_nmds, smp_traits_scaled, permutations = 999)
fame_trait_correlations <- envfit(fame_nmds, fame_traits_scaled, permutations = 999)

# Extract significant vectors
smp_significant_vectors <- as.data.frame(smp_trait_correlations$vectors$arrows)
smp_scaled_vectors <- smp_significant_vectors %>%
  mutate(NMDS1 = NMDS1 * 0.5,
         NMDS2 = NMDS2 * 0.5)
rownames(smp_significant_vectors) <- rownames(smp_trait_correlations$vectors$arrows)

fame_significant_vectors <- as.data.frame(fame_trait_correlations$vectors$arrows)
fame_scaled_vectors <- fame_significant_vectors %>%
  mutate(NMDS1 = NMDS1 * 0.5,
         NMDS2 = NMDS2 * 0.5)
rownames(fame_significant_vectors) <- rownames(fame_trait_correlations$vectors$arrows)

# Generate NMDS plot with trait vectors
fame_nmds_scores$NMDS1 <- fame_nmds_scores[,1]
fame_nmds_scores$NMDS2 <- fame_nmds_scores[,2]
fame.mammal.pcoa.vectors.plot <- plot_pcoa_with_vectors(fame_nmds_scores, fame_scaled_vectors)
smp_nmds_scores$NMDS1 <- smp_nmds_scores[,1]
smp_nmds_scores$NMDS2 <- smp_nmds_scores[,2]
smp.mammal.pcoa.vectors.plot <- plot_pcoa_with_vectors(smp_nmds_scores, smp_scaled_vectors)

(mammal.pcoa.vectors.plot = plot_grid(fame.mammal.pcoa.vectors.plot, smp.mammal.pcoa.vectors.plot, labels=c("a)","b)")))

ggsave(plot = mammal.pcoa.vectors.plot , filename = "plots/mammal.pcoa.vectors.plot.jpg", width = 8, height = 3.5, scale =1.5)

# Save
write.csv(mammal_traits, "data_clean/mammal_traits_clustered.csv")

#################
#### Birds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
#### Find Important Traits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Part 1 - Select the traits for inclusion
birds_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Birds") %>% mutate(stratum = as.factor(stratum),
                                                                                   stratum_aquatic = as.factor(stratum_aquatic),
                                                                                      nesting = as.factor(nesting),
                                                                                      diet= as.factor(diet),
                                                                                      dominant_pyrome = as.factor(dominant_pyrome)) %>% droplevels()
# Traits of interest
all_traits <- c("Mass_g", 
                "stratum_aerial", "stratum_arboreal_insessorial", "stratum_aquatic", "stratum_terrestrial", 
                "nest_ground", "nest_hollows", "nest_branch", 
                "diet_invertivore", "diet_omnivore", 
                "litter_size_n", "litters_per_year_n", "max_longevity_d", "n_offspring_year",
                "biggest_patch_size", "n_habitat_patches", 
                "dominant_pyrome", "pyrome_breadth")
# Create the formula
response_var <- "meanben_firefreq"
traits_formula <- as.formula(paste(response_var, "~", paste(all_traits, collapse = " + ")))

#### Too Frequent Fire
# Fit the model as a Boosted Regression Tree
smp_modelling = birds_traits %>% select(meanben_firefreq, all_of(all_traits)) %>%
  drop_na(meanben_firefreq)

birds.brt.smp <- gbm.step(data = smp_modelling, 
                          gbm.x = 2:19,
                          gbm.y = 1,
                          family = "gaussian", 
                          tree.complexity = 2, 
                          learning.rate = 0.001,
                          bag.fraction = 0.5)
# Look at the variable importance
smp.bird.importance = summary(birds.brt.smp)  # Variable importance

smp.bird.importance.plot = smp.bird.importance %>%
  left_join(trait_df, by = "var") %>%
  mutate(traits_label = fct_reorder(traits_label, rel.inf, .desc = FALSE)) %>%  # Reorder `var` based on `rel.inf`
  mutate(fire.var = "Frequent Fire")
  
# Partial Dependence Plots
factor_vars <- names(smp_modelling)[sapply(smp_modelling, is.factor)==TRUE]
pdp_data <- lapply(birds.brt.smp$var.names, get_pdp_data, model =birds.brt.smp)
bird_pdps <- lapply(pdp_data, create_pdp_plots, factor_vars)
bird_smp_pdp <- wrap_plots(bird_pdps, ncol = 4)  # 5 rows, 3 columns
bird_smp_pdp
ggsave(plot = bird_smp_pdp, filename = "plots/bird.smp.pdp.jpg", width = 8, height = 5, scale=2)

#### Recent Fire
fame_modelling = birds_traits %>% select(fame_lm_slope, all_of(all_traits)) %>%
  drop_na(fame_lm_slope)%>%
  mutate(fame_lm_slope = 1000*fame_lm_slope)

birds.brt.fame <- gbm.step(data = fame_modelling, 
                          gbm.x = 2:19,
                          gbm.y = 1,
                          family = "gaussian", 
                          tree.complexity = 2, 
                          learning.rate = 0.001,
                          bag.fraction = 0.5)

# Partial Dependence Plots
factor_vars <- names(fame_modelling)[sapply(fame_modelling, is.factor)==TRUE]
pdp_data <- lapply(birds.brt.fame$var.names, get_pdp_data, model =birds.brt.fame)
bird_pdps <- lapply(pdp_data, create_pdp_plots, factor_vars)
bird_fame_pdp <- wrap_plots(bird_pdps, ncol = 4)  # 5 rows, 3 columns
bird_fame_pdp
ggsave(plot = bird_fame_pdp, filename = "plots/bird.fame.pdp.jpg", width = 8, height = 5, scale=2)

# Look at the variable importance
fame.bird.importance = summary(birds.brt.fame)  # Variable importance
fame.bird.importance.plot = fame.bird.importance %>%
  left_join(trait_df, by = "var") %>%
  mutate(traits_label = fct_reorder(traits_label, rel.inf, .desc = FALSE)) %>%  # Reorder `var` based on `rel.inf`
  mutate(fire.var = "Recent Fire")

# Plot importance
bird.importance = rbind(smp.bird.importance.plot, fame.bird.importance.plot)

bird.importance.plot = bird.importance %>%
  filter(rel.inf > 0) %>%
  mutate(fire.var = factor(fire.var, levels = c("Recent Fire", "Frequent Fire"))) %>%
ggplot(aes(x = traits_label, y = rel.inf, fill = Group)) + geom_col() + 
  labs(fill = "Trait Group", x="Trait", y = "Relative Influence (%)", title="Bird Trait Importance") + 
  facet_wrap(~fire.var, nrow=1) +
  scale_fill_viridis_d() + theme_minimal() +
  coord_flip()

bird.importance.plot
ggsave(plot = bird.importance.plot, filename = "plots/combined.bird.importance.plot.jpg", width = 8, height = 3.5)

#### Clustering Species by Important Traits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
smp_traits_scaled = birds_traits %>% 
  select(all_of(smp.bird.importance$var))  %>% 
  select(all_of(smp.bird.importance$var)) # Order to match with vector of importance
smp.bird.trait.importance = smp.bird.importance$rel.inf/100

# Compute Gower distance matrix
smp_gower_dist <- daisy(smp_traits_scaled, metric = "gower", weights=smp.bird.trait.importance)

smp_gmm_model <- Mclust(smp_gower_dist, G=1:(nrow(smp_traits_scaled)/5)) # Try between 1-20 clusters to see which is best. 
summary(smp_gmm_model)

# Add cluster assignments to the dataset
birds_traits$smp_cluster <- smp_gmm_model$classification
# Summarize meanben_firefreq by cluster
species_clusters <- birds_traits %>% select(Common_Name, smp_cluster)

fame_traits_scaled = birds_traits %>% 
  select(all_of(fame.bird.importance$var))  %>% 
  select(all_of(fame.bird.importance$var)) # Order to match with vector of importance
fame.bird.trait.importance = fame.bird.importance$rel.inf/100

# Compute Gower distance matrix
fame_gower_dist <- daisy(fame_traits_scaled, metric = "gower", weights=fame.bird.trait.importance)

fame_gmm_model <- Mclust(fame_gower_dist, G=1:(nrow(fame_traits_scaled)/5)) # Try between 1-20 clusters to see which is best. 
summary(fame_gmm_model)

# Add cluster assignments to the dataset
birds_traits$fame_cluster <- fame_gmm_model$classification
# Summarize meanben_firefreq by cluster
species_clusters <- birds_traits %>% select(Common_Name, fame_cluster)

## NMDS
smp_nmds <- metaMDS(smp_gower_dist,method="gower", k = 2, trymax = 30)
# Extract NMDS scores and perform clustering
smp_nmds_scores <- as.data.frame(scores(smp_nmds))
smp_nmds_scores$Species = birds_traits$Common_Name
birds_traits$smp_cluster = as.factor(birds_traits$smp_cluster)
smp_nmds_scores$Cluster = as.factor(birds_traits$smp_cluster)
smp_nmds_scores$aquatic = as.factor(birds_traits$stratum_aquatic)

# Recent fire
fame_nmds <- metaMDS(fame_gower_dist,method="gower", k = 2, trymax = 30)
# Extract NMDS scores and perform clustering
fame_nmds_scores <- as.data.frame(scores(fame_nmds))
fame_nmds_scores$Species = birds_traits$Common_Name
birds_traits$fame_cluster = as.factor(birds_traits$fame_cluster)
fame_nmds_scores$Cluster = as.factor(birds_traits$fame_cluster)
fame_nmds_scores$aquatic = as.factor(birds_traits$stratum_aquatic)

# Plot NMDS with clusters
(smp.bird.pcoa = plot_pcoa_aquatic(smp_nmds_scores))
(fame.bird.pcoa = plot_pcoa_aquatic(fame_nmds_scores))

combined.bird.pcoa <- plot_grid(fame.bird.pcoa, smp.bird.pcoa, nrow=2)

ggsave(plot = combined.bird.pcoa, filename = "plots/combined.bird.pcoa.plot.jpg", width = 8, height = 7, scale =2)

write.csv(birds_traits, "data_clean/birds_traits_clustered.csv")

#################
#### Reptiles #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
reptile_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Reptiles") %>%
    mutate(across(c(stratum, nesting, diet, dominant_pyrome, stratum_aerial, stratum_arboreal_insessorial, stratum_aquatic,
                    stratum_fossorial, stratum_saxicolous, stratum_terrestrial), as.factor))
# Traits of interest
all_traits <- c("Mass_g", "home_range_km2", "dispersal_km", 
                "stratum_arboreal_insessorial", "stratum_aquatic",
                "stratum_fossorial", "stratum_saxicolous", "stratum_terrestrial",
                "diet_carnivore", "diet_invertivore", "diet_omnivore",
                "litter_size_n", "litters_per_year_n", "max_longevity_d", "n_offspring_year",
                "biggest_patch_size", "n_habitat_patches", "dominant_pyrome", "pyrome_breadth")

# Create the formula
response_var <- "meanben_firefreq"
smp_traits_formula <- as.formula(paste(response_var, "~", paste(all_traits, collapse = " + ")))
response_var <- "fame_lm_slope"
fame_traits_formula <- as.formula(paste(response_var, "~", paste(all_traits, collapse = " + ")))

# Fit the model as a Boosted Regression Tree
smp_modelling = reptile_traits  %>% select(meanben_firefreq, all_of(all_traits)) %>%
  drop_na(meanben_firefreq)
fame_modelling = reptile_traits  %>% select(fame_lm_slope, all_of(all_traits)) %>%
  drop_na(fame_lm_slope)

# Recent Fire
reptile.brt.smp <- gbm(
  formula = smp_traits_formula,
  data = smp_modelling,
  distribution = "gaussian",
  n.trees = 3000,        # Number of trees
  interaction.depth = 1, # Tree depth to capture interactions
  shrinkage = 0.00001,      # Learning rate
  bag.fraction = 0.5, 
  cv.folds = 10           # Cross-validation for tuning
)

(smp.reptile.importance <- summary(reptile.brt.smp))

reptile.brt.fame <- gbm(
  formula = fame_traits_formula,
  data = fame_modelling,
  distribution = "gaussian",
  n.trees = 3000,        # Number of trees
  interaction.depth = 1, # Tree depth to capture interactions
  shrinkage = 0.00001,      # Learning rate
  bag.fraction = 0.5, 
  cv.folds = 10           # Cross-validation for tuning
)

(fame.reptile.importance <- summary(reptile.brt.fame))


# Partial Dependence Plots
factor_vars <- names(fame_modelling)[sapply(fame_modelling, is.factor)==TRUE]
pdp_data <- lapply(reptile.brt.fame$var.names, get_pdp_data, model =reptile.brt.fame)
reptile_pdps <- lapply(pdp_data, create_pdp_plots, factor_vars)
reptile_fame_pdp <- wrap_plots(reptile_pdps, ncol = 4)  # 5 rows, 3 columns
reptile_fame_pdp
ggsave(plot = reptile_fame_pdp, filename = "plots/reptile.fame.pdp.jpg", width = 8, height = 5, scale=2)

factor_vars <- names(smp_modelling)[sapply(smp_modelling, is.factor)==TRUE]
pdp_data <- lapply(reptile.brt.smp$var.names, get_pdp_data, model =reptile.brt.smp)
reptile_pdps <- lapply(pdp_data, create_pdp_plots, factor_vars)
reptile_smp_pdp <- wrap_plots(reptile_pdps, ncol = 4)  # 5 rows, 3 columns
reptile_smp_pdp
ggsave(plot = reptile_smp_pdp, filename = "plots/reptile.smp.pdp.jpg", width = 8, height = 5, scale=2)

# Look at the variable importance
fame.reptile.importance.plot = fame.reptile.importance %>%
  left_join(trait_df, by = "var") %>%
  mutate(var = fct_reorder(var, rel.inf, .desc = FALSE)) %>%  # Reorder `var` based on `rel.inf`
  mutate(fire.var = "Recent Fire")

smp.reptile.importance.plot = smp.reptile.importance %>%
  left_join(trait_df, by = "var") %>%
  mutate(var = fct_reorder(var, rel.inf, .desc = FALSE)) %>%  # Reorder `var` based on `rel.inf`
  mutate(fire.var = "Frequent Fire")

reptile.importance = rbind(smp.reptile.importance.plot, fame.reptile.importance.plot)

reptile.importance.plot = reptile.importance %>%   
  filter(rel.inf > 0) %>%
  mutate(fire.var = factor(fire.var, levels = c("Recent Fire", "Frequent Fire"))) %>%
ggplot(aes(x = var, y = rel.inf, fill = Group)) + geom_col() + 
  labs(fill = "Trait Group", x="Trait", y = "Relative Influence (%)", title="Reptile Trait Importance") + 
  scale_fill_viridis_d() + theme_minimal() + facet_wrap(~fire.var,nrow=1) +
  coord_flip()

reptile.importance.plot
ggsave(plot = reptile.importance.plot, filename = "plots/combined.reptile.importance.plot.jpg", width = 8, height = 4)

####################################################
#### Clustering Species by Important Traits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################
### Recent fire
fame_traits_scaled = reptile_traits %>% 
  select(all_of(fame.reptile.importance$var))  %>% 
  select(all_of(fame.reptile.importance$var)) # Order to match with vector of importance

fame.trait.importance = fame.reptile.importance$rel.inf/100

# Compute Gower distance matrix
fame_gower_dist <- daisy(fame_traits_scaled, metric = "gower", 
                    weights=fame.trait.importance)

fame_gmm_model <- Mclust(fame_gower_dist, G=1:22) # Try between 1-20 clusters to see which is best. 

# Add cluster assignments to the dataset
reptile_traits$fame_cluster <- fame_gmm_model$classification

# Summarize meanben_firefreq by cluster
species_clusters <- reptile_traits %>% select(Common_Name, fame_cluster)

## NMDS
fame_nmds <- metaMDS(fame_gower_dist,method="gower", k = 2, trymax = 30)

# Extract NMDS scores and perform clustering
fame_nmds_scores <- as.data.frame(scores(fame_nmds))
fame_nmds_scores$Species = reptile_traits$Common_Name
fame_nmds_scores$Cluster = as.factor(reptile_traits$fame_cluster)

### Frequent fire
smp_traits_scaled = reptile_traits %>% 
  select(all_of(smp.reptile.importance$var))  %>% 
  select(all_of(smp.reptile.importance$var)) # Order to match with vector of importance

smp.trait.importance = smp.reptile.importance$rel.inf/100

# Compute Gower distance matrix
smp_gower_dist <- daisy(smp_traits_scaled, metric = "gower", 
                         weights=smp.trait.importance)

smp_gmm_model <- Mclust(smp_gower_dist, G=1:22) # Try between 1-20 clusters to see which is best. 

# Add cluster assignments to the dataset
reptile_traits$smp_cluster <- smp_gmm_model$classification

# Summarize meanben_firefreq by cluster
species_clusters <- reptile_traits %>% select(Common_Name, smp_cluster)

## NMDS
smp_nmds <- metaMDS(smp_gower_dist,method="gower", k = 2, trymax = 30)

# Extract NMDS scores and perform clustering
smp_nmds_scores <- as.data.frame(scores(smp_nmds))
smp_nmds_scores$Species = reptile_traits$Common_Name
smp_nmds_scores$Cluster = as.factor(reptile_traits$smp_cluster)

# Plot NMDS with clusters
(fame.reptile.pcoa = plot_pcoa(fame_nmds_scores))
(smp.reptile.pcoa = plot_pcoa(smp_nmds_scores))

(combined.reptile.pcoa = plot_grid(fame.reptile.pcoa, smp.reptile.pcoa, labels=c("a)", "b)")))

ggsave(plot = combined.reptile.pcoa, filename = "plots/combined.reptile.pcoa.plot.jpg", width = 8, height = 3.5, scale =2)

write.csv(reptile_traits, "data_clean/reptile_traits_clustered.csv")

#################
#### Frogs #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
frog_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Amphibians") %>% mutate(stratum = as.factor(stratum),
                                                                                       stratum_aquatic = as.factor(stratum_aquatic),
                                                                                       stratum_arboreal_insessorial = as.factor(stratum_arboreal_insessorial),
                                                                                       stratum_fossorial = as.factor(stratum_fossorial),
                                                                                       nest_burrow = as.factor(nest_burrow),
                                                                                      nesting = as.factor(nesting),
                                                                                      diet= as.factor(diet),
                                                                                      dominant_pyrome = as.factor(dominant_pyrome))

# Traits of interest
all_traits <- c("scale(Mass_g)",
                "scale(litter_size_n)",
                "stratum_arboreal_insessorial", "stratum_aquatic",  "stratum_fossorial", "nest_burrow", 
                "scale(biggest_patch_size)", "scale(n_habitat_patches)", "dominant_pyrome", "scale(pyrome_breadth)")
# Create the formula
response_var <- "meanben_firefreq"
smp_formula <- as.formula(paste(response_var, "~", paste(all_traits, collapse = " + ")))
dat_modelling = frog_traits %>% drop_na(meanben_firefreq)
frog.model = glm(smp_formula, data = dat_modelling)

# Frequent Fire
factor_vars <- names(dat_modelling)[sapply(dat_modelling, is.factor)==TRUE]
pred_data <- lapply(str_remove_all(all_traits, "scale\\(|\\)"), get_pred_data, model =frog.model)
frog_preds <- lapply(pred_data, create_pred_plots, factor_vars)
frog_smp_preds <- wrap_plots(frog_preds, ncol = 4)  # 5 rows, 3 columns
frog_smp_preds
ggsave(plot = reptile_smp_pdp, filename = "plots/reptile.smp.pdp.jpg", width = 8, height = 5, scale=2)


brt_model <- gbm(
  formula = smp_formula,
  data = dat_modelling,
  distribution = "gaussian",
  n.trees = 3000,        # Number of trees
  interaction.depth = 1, # Tree depth to capture interactions
  shrinkage = 0.01,      # Learning rate
  bag.fraction = 0.5, 
  cv.folds = 10           # Cross-validation for tuning
)


####################################################
#### Clustering Species by Important Traits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################
frog_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Amphibians") %>% mutate(stratum = as.factor(stratum),
                                                                                       nesting = as.factor(nesting),
                                                                                       diet= as.factor(diet),
                                                                                       stratum_aquatic = as.factor(stratum_aquatic),
                                                                                       stratum_arboreal_insessorial = as.factor(stratum_arboreal_insessorial),
                                                                                       stratum_fossorial = as.factor(stratum_fossorial),
                                                                                       nest_burrow = as.factor(nest_burrow),
                                                                                       dominant_pyrome = as.factor(dominant_pyrome))
all_traits <- c("Mass_g","stratum_arboreal_insessorial", "stratum_aquatic", "stratum_aquatic",  "stratum_fossorial", "nest_burrow", "litter_size_n", "biggest_patch_size", "n_habitat_patches", "dominant_pyrome", "pyrome_breadth")

traits_scaled = frog_traits %>% select(Common_Name,meanben_firefreq, all_of(all_traits)) %>% drop_na()

# Compute Gower distance matrix
gower_dist <- daisy(traits_scaled[,3:12], metric = "gower")

gmm_model <- Mclust(gower_dist, G=1:7) # Try between 1-20 clusters to see which is best. 

# Summarize the model
summary(gmm_model)

# Plot the clustering result
plot(gmm_model, what = "BIC")
#plot(gmm_model, what = "classification")

# Add cluster assignments to the dataset
traits_scaled$cluster <- gmm_model$classification

frog_traits = left_join(frog_traits, traits_scaled)

# Summarize meanben_firefreq by cluster
species_clusters <- traits_scaled %>% select(Common_Name, cluster)

trait_summary <- traits_scaled %>%
  group_by(cluster) %>%
  summarise(mean_firefreq = mean(meanben_firefreq, na.rm = TRUE),
            sd_firefreq = sd(meanben_firefreq, na.rm = TRUE))

ggplot(traits_scaled) + geom_boxplot(aes(x=as.factor(cluster), y = meanben_firefreq))

## NMDS
nmds <- metaMDS(gower_dist,method="gower", k = 2, trymax = 30)
# Extract NMDS scores and perform clustering
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Species = traits_scaled$Common_Name
nmds_scores$Cluster = as.factor(traits_scaled$cluster)

# Plot NMDS with clusters
(frog.pcoa <- plot_pcoa(nmds_scores))

ggsave(plot = frog.pcoa, filename = "plots/frog.pcoa.plot.jpg", width = 6, height = 4, scale =2)

write.csv(frog_traits, "data_clean/frog_traits_clustered.csv")

#### ENDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

