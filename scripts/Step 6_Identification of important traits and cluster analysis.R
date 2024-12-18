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
#          with an extra column specifying the cluster in which each species has been assigned
#          based on the cluster analysis. 
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

source("scripts/brt.functions.R")

# Read in the imputed data
vic_fauna_traits_imputed = readRDS("data_clean/vic_fauna_traits_imputed.Rds")

#################
#### Mammals #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
#### Find Important Traits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mammal_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Mammals") %>% mutate(stratum = as.factor(stratum),
                                                                                      nesting = as.factor(nesting),
                                                                                      diet = as.factor(diet),
                                                                                      dominant_pyrome = as.factor(dominant_pyrome),
                                                                                      volant = as.factor(volant), 
                                                                                      hibernation_torpor = as.factor(hibernation_torpor))

# Define the trait groups
trait_groups = list(Movement_Physiology = c("Mass_g", "home_range_km2", "dispersal_km",  "max_longevity_d"),
                    Reproduction = c("litter_size_n", "litters_per_year_n","n_offspring_year"),
                    Broad_Habitat = c("dominant_pyrome"),
                    Microhabitat = c("stratum","stratum_aerial","stratum_arboreal_insessorial", "stratum_aquatic","stratum_cryptic","stratum_fossorial", "stratum_terrestrial", "stratum_saxicolous"),
                    Nesting = c("nesting", "nest_burrow","nest_cave","nest_ground","nest_hollows","nest_branch"),
                    Behaviour = c("hibernation_torpor", "volant"),
                    Diet = c("diet","diet_breadth_n","diet_animals","diet_plants", "diet_carnivore","diet_herbivore",
                             "diet_invertivore","diet_infloresence","diet_omnivore","diet_granivore", "diet_simple"),
                    Biogeography = c("biggest_patch_size","n_habitat_patches","patch_isolation", "pyrome_breadth"))

# Convert the list to a data frame using stack
trait_df <- stack(trait_groups)

# Rename columns for clarity
colnames(trait_df) <- c("var", "Group")

# Traits of interest
all_traits <- c("Mass_g", "home_range_km2", "dispersal_km", "stratum", "nesting", "volant", "hibernation_torpor", 
                "diet", "litter_size_n", "litters_per_year_n", "max_longevity_d", "n_offspring_year",
                "biggest_patch_size", "n_habitat_patches", "dominant_pyrome", "pyrome_breadth")

all_traits_clean_names = data.frame(var = all_traits, 
                                    traits_label = c("Mass (g)", "Home Range (km2)", "Dispersal (km)", "Stratum", "Nesting Preference",
                                                     "Volant", "Hibernation/Torpor", "Diet Category", "Litter Size", "Litters/Yr", "Longevity", "Offspring/Yr", 
                                                     "Size of Largest Habitat Patch", "Number of Habitat Patches", "Dominant Pyrome", "Pyrome Breadth"))
trait_df = left_join(trait_df, all_traits_clean_names)

# Create the formula
traits_formula_smp <- as.formula(paste("meanben_firefreq", "~", paste(all_traits, collapse = " + ")))
traits_formula_fame <- as.formula(paste("fame_lm_slope", "~", paste(all_traits, collapse = " + ")))

# Fit the model as a Boosted Regression Tree
smp_modelling = mammal_traits %>% select(meanben_firefreq, all_of(all_traits)) %>%
  drop_na(meanben_firefreq)

fame_modelling = mammal_traits %>% select(fame_lm_slope, all_of(all_traits)) %>%
  drop_na(fame_lm_slope)

mammal_brt_smp <- gbm.step(data = smp_modelling, 
                         gbm.x = 2:17,
                         gbm.y = 1,
                         family = "gaussian", 
                         tree.complexity = 2, 
                         learning.rate = 0.005,
                         bag.fraction = 0.5)

summary(mammal_brt_smp)
plot.gbm(mammal_brt_smp, i.var = 8)

# Look at the variable importance
mammal.importance = summary(mammal_brt_smp)  # Variable importance
all_traits_clean_names$traits_label = as.factor(all_traits_clean_names$traits_label)

mammal.importance.plot = mammal.importance %>%
  left_join(trait_df, by = "var") %>%
  mutate(traits_label = fct_reorder(traits_label, rel.inf, .desc = FALSE)) %>%  # Reorder `var` based on `rel.inf`
  ggplot(aes(x = traits_label, y = rel.inf, fill = Group)) + geom_col() + 
  labs(fill = "Trait Group", x="Trait", y = "Relative Influence (%)", title="Mammal Trait Importance") + 
  scale_fill_viridis_d() + theme_minimal() +
  coord_flip()

ggsave(plot = mammal.importance.plot, filename = "plots/mammal.importance.plot.jpg", width = 8, height = 4)

#### Clustering Species by Important Traits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
traits_scaled = mammal_traits %>% 
  select(all_of(all_traits)) %>% 
  select(all_of(mammal.importance$var)) # Order to match with vector of importance

mammal.importance = mammal.importance$rel.inf/100

# Compute Gower distance matrix
gower_dist <- daisy(traits_scaled, metric = "gower", weights = mammal.importance)
gmm_model <- Mclust(gower_dist, G=1:12) # Try between 1-20 clusters to see which is best. 

# Add cluster assignments to the dataset
mammal_traits$cluster <- gmm_model$classification

# Summarize meanben_firefreq by cluster
species_clusters <- mammal_traits %>% select(Common_Name, cluster)

trait_summary <- mammal_traits %>%
  group_by(cluster) %>%
  summarise(mean_firefreq = mean(meanben_firefreq, na.rm = TRUE),
            sd_firefreq = sd(meanben_firefreq, na.rm = TRUE))

ggplot(mammal_traits) + geom_boxplot(aes(x=as.factor(cluster), y = meanben_firefreq))

## NMDS
nmds <- metaMDS(gower_dist,method="gower", k = 2, trymax = 30)

# Extract NMDS scores and perform clustering
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Species = mammal_traits$Common_Name
nmds_scores$Cluster = as.factor(mammal_traits$cluster)

# Plot NMDS with clusters
plot_pcoa <- function(data.in) {
  # Calculate convex hulls for each group
  hulls <- data.in %>%
    group_by(Cluster) %>%  # Use tidy evaluation for dynamic grouping
    slice(chull(NMDS1, NMDS2))           # Select points forming the convex hull
  
  # Create the PCA plot
  pcoa_plot <- ggplot(data.in) + 
    geom_point(aes(x = NMDS1, y = NMDS2, colour = Cluster), size = 4) + 
    geom_text(aes(x = NMDS1, y = NMDS2, label = Species), nudge_y = -0.01, size = 3) +
    geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, 
                                   group = Cluster, 
                                   fill = Cluster), alpha = 0.3) +
    theme_bw() +
    xlab("NMDS1") + 
    ylab("NMDS2") #+ facet_wrap(~Cluster, scales="free")
  
  return(pcoa_plot)
}

(mammal.pcoa = plot_pcoa(nmds_scores))
ggsave(plot = mammal.pcoa, filename = "plots/mammal.pcoa.plot.jpg", width = 6, height = 4, scale =2)


## Add the traits to the plot
# Compute correlations of traits with NMDS axes
trait_correlations <- envfit(nmds, traits_scaled, permutations = 999)

# Extract significant vectors
significant_vectors <- as.data.frame(trait_correlations$vectors$arrows)
scaled_vectors <- significant_vectors %>%
  mutate(NMDS1 = NMDS1 * 0.5,
         NMDS2 = NMDS2 * 0.5)
rownames(significant_vectors) <- rownames(trait_correlations$vectors$arrows)

# Add to NMDS plot
plot_pcoa_with_vectors <- function(data.in, vectors) {
  # Calculate convex hulls for each group
  hulls <- data.in %>%
    group_by(Cluster) %>% 
    slice(chull(NMDS1, NMDS2))
  
  # Create NMDS plot with vectors
  pcoa_plot <- ggplot(data.in) + 
    geom_point(aes(x = NMDS1, y = NMDS2, colour = Cluster), size = 4) + 
    geom_segment(data = vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
    geom_text(data = vectors, aes(x = NMDS1, y = NMDS2, label = rownames(vectors)),
              color = "black", hjust = 0, vjust = -0.5, position=position_dodge2(width=1)) + 
    coord_cartesian(ylim=c(-0.6,0.6), xlim=c(-0.5, 0.9)) +
    theme_bw() +
    xlab("NMDS1") + 
    ylab("NMDS2")
  
  return(pcoa_plot)
}

# Generate NMDS plot with trait vectors
nmds_scores$NMDS1 <- nmds_scores[,1]
nmds_scores$NMDS2 <- nmds_scores[,2]
mammal.pcoa.vectors.plot <- plot_pcoa_with_vectors(nmds_scores, scaled_vectors)
ggsave(plot = mammal.pcoa.vectors.plot , filename = "plots/mammal.pcoa.vectors.plot.jpg", width = 6, height = 4, scale =2)

# Save
write.csv(mammal_traits, "data_clean/mammal_traits_clustered.csv")

#################
#### Birds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
#### Find Important Traits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
birds_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Birds") %>% mutate(stratum = as.factor(stratum),
                                                                                   stratum_aquatic = as.factor(stratum_aquatic),
                                                                                      nesting = as.factor(nesting),
                                                                                      diet= as.factor(diet),
                                                                                      dominant_pyrome = as.factor(dominant_pyrome)) %>% droplevels()
# Traits of interest
all_traits <- c("Mass_g", "stratum", "nesting", "stratum_aquatic",
                "diet", "litter_size_n", "litters_per_year_n", "max_longevity_d", "n_offspring_year",
                "biggest_patch_size", "n_habitat_patches", 
                "dominant_pyrome", "pyrome_breadth")
# Create the formula
response_var <- "meanben_firefreq"
traits_formula <- as.formula(paste(response_var, "~", paste(all_traits, collapse = " + ")))

# Fit the model as a Boosted Regression Tree
smp_modelling = birds_traits %>% select(meanben_firefreq, all_of(all_traits)) %>%
  drop_na(meanben_firefreq)

birds.brt.smp <- gbm.step(data = smp_modelling, 
                         gbm.x = 2:14,
                         gbm.y = 1,
                         family = "gaussian", 
                         tree.complexity = 2, 
                         learning.rate = 0.005,
                         bag.fraction = 0.5)

summary(birds.brt.smp)
plot.gbm(birds.brt.smp, i.var = 5)

# Look at the variable importance
bird.importance = summary(birds.brt.smp)  # Variable importance

bird.importance.plot = bird.importance %>%
  left_join(trait_df, by = "var") %>%
  mutate(var = fct_reorder(var, rel.inf, .desc = FALSE)) %>%  # Reorder `var` based on `rel.inf`
  ggplot(aes(x = var, y = rel.inf, fill = Group)) + geom_col() + 
  labs(fill = "Trait Group", x="Trait", y = "Relative Influence (%)", title="Mammal Trait Importance") + 
  scale_fill_viridis_d() + theme_minimal() +
  coord_flip()

ggsave(plot = bird.importance.plot, filename = "plots/bird.importance.plot.jpg", width = 8, height = 4)

#### Clustering Species by Important Traits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
traits_scaled = birds_traits %>% 
  select(all_of(bird.importance$var))  %>% 
  select(all_of(bird.importance$var)) # Order to match with vector of importance

bird.trait.importance = bird.importance$rel.inf/100

# Compute Gower distance matrix
gower_dist <- daisy(traits_scaled, metric = "gower", weights=bird.trait.importance)

gmm_model <- Mclust(gower_dist, G=1:25) # Try between 1-20 clusters to see which is best. 

# Summarize the model
summary(gmm_model)

# Add cluster assignments to the dataset
birds_traits$cluster <- gmm_model$classification

# Summarize meanben_firefreq by cluster
species_clusters <- birds_traits %>% select(Common_Name, cluster)

trait_summary <- birds_traits %>%
  group_by(cluster) %>%
  summarise(mean_firefreq = mean(meanben_firefreq, na.rm = TRUE),
            sd_firefreq = sd(meanben_firefreq, na.rm = TRUE))
ggplot(birds_traits) + geom_boxplot(aes(x=as.factor(cluster), y = meanben_firefreq))

## NMDS
nmds <- metaMDS(gower_dist,method="gower", k = 2, trymax = 30)
# Extract NMDS scores and perform clustering
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Species = birds_traits$Common_Name
birds_traits$cluster = as.factor(birds_traits$cluster)
nmds_scores$Cluster = as.factor(birds_traits$cluster)
nmds_scores$aquatic = as.factor(birds_traits$stratum_aquatic)

# Plot NMDS with clusters
plot_pcoa <- function(data.in) {
  # Calculate convex hulls for each group
  hulls <- data.in %>%
    group_by(Cluster) %>%  # Use tidy evaluation for dynamic grouping
    slice(chull(NMDS1, NMDS2))           # Select points forming the convex hull
  
  # Create the PCA plot
  pcoa_plot <- ggplot(data.in) + 
    geom_point(aes(x = NMDS1, y = NMDS2, colour = Cluster), size = 4) + 
    geom_text(aes(x = NMDS1, y = NMDS2, label = Species), nudge_y = -0.01) +
    geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, 
                                   group = Cluster, 
                                   fill = Cluster), alpha = 0.3) +
    theme_bw() + #facet_wrap(~Cluster, scales="free", nrow=2) +
    xlab("NMDS1") + 
    ylab("NMDS2") + facet_wrap(~aquatic, scales = "free")
  
  return(pcoa_plot)
}

(bird.pcoa = plot_pcoa(nmds_scores))
ggsave(plot = bird.pcoa, filename = "plots/bird.pcoa.plot.jpg", width = 8, height = 4, scale =2)

write.csv(birds_traits, "data_clean/birds_traits_clustered.csv")

#################
#### Reptiles #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
reptile_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Reptiles") %>%
    mutate(across(c(stratum, nesting, diet, dominant_pyrome, stratum_aerial, stratum_arboreal_insessorial, stratum_aquatic,
                    stratum_fossorial, stratum_saxicolous, stratum_terrestrial), as.factor))
# Traits of interest
all_traits <- c("Mass_g", "home_range_km2", "dispersal_km", "stratum_arboreal_insessorial", "stratum_aquatic",
                "stratum_fossorial", "stratum_saxicolous", "stratum_terrestrial",
                "diet", "litter_size_n", "litters_per_year_n", "max_longevity_d", "n_offspring_year",
                "biggest_patch_size", "n_habitat_patches", "dominant_pyrome", "pyrome_breadth")

# Create the formula
response_var <- "meanben_firefreq"
traits_formula <- as.formula(paste(response_var, "~", paste(all_traits, collapse = " + ")))

# Fit the model as a Boosted Regression Tree
smp_modelling = reptile_traits  %>% select(meanben_firefreq, all_of(all_traits)) %>%
  drop_na(meanben_firefreq)

reptile.brt.smp <- gbm(
  formula = traits_formula,
  data = smp_modelling,
  distribution = "gaussian",
  n.trees = 3000,        # Number of trees
  interaction.depth = 1, # Tree depth to capture interactions
  shrinkage = 0.01,      # Learning rate
  bag.fraction = 0.5, 
  cv.folds = 10           # Cross-validation for tuning
)

(reptile.importance <- summary(reptile.brt.smp))

# Look at the variable importance
reptile.importance = summary(reptile.brt.smp)  # Variable importance
reptile.importance.plot = reptile.importance %>%
  left_join(trait_df, by = "var") %>%
  mutate(var = fct_reorder(var, rel.inf, .desc = FALSE)) %>%  # Reorder `var` based on `rel.inf`
  ggplot(aes(x = var, y = rel.inf, fill = Group)) + geom_col() + 
  labs(fill = "Trait Group", x="Trait", y = "Relative Influence (%)", title="Reptile Trait Importance") + 
  scale_fill_viridis_d() + theme_minimal() +
  coord_flip()

ggsave(plot = reptile.importance.plot, filename = "plots/reptile.importance.plot.jpg", width = 8, height = 4)

####################################################
#### Clustering Species by Important Traits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################################################

traits_scaled = reptile_traits %>% 
  select(all_of(reptile.importance$var))  %>% 
  select(all_of(reptile.importance$var)) # Order to match with vector of importance

trait.importance = reptile.importance$rel.inf/100

# Compute Gower distance matrix
gower_dist <- daisy(traits_scaled, metric = "gower", 
                    #weights=trait.importance
                    )

gmm_model <- Mclust(gower_dist, G=1:25) # Try between 1-20 clusters to see which is best. 

# Summarize the model
summary(gmm_model)

# Plot the clustering result
plot(gmm_model, what = "BIC")
#plot(gmm_model, what = "classification")

# Add cluster assignments to the dataset
reptile_traits$cluster <- gmm_model$classification

# Summarize meanben_firefreq by cluster
species_clusters <- reptile_traits %>% select(Common_Name, cluster)

trait_summary <- reptile_traits %>%
  group_by(cluster) %>%
  summarise(mean_firefreq = mean(meanben_firefreq, na.rm = TRUE),
            sd_firefreq = sd(meanben_firefreq, na.rm = TRUE))

ggplot(reptile_traits) + geom_boxplot(aes(x=as.factor(cluster), y = meanben_firefreq))

## NMDS
nmds <- metaMDS(gower_dist,method="gower", k = 2, trymax = 30)

# Extract NMDS scores and perform clustering
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Species = reptile_traits$Common_Name
nmds_scores$Cluster = as.factor(reptile_traits$cluster)

# Plot NMDS with clusters
plot_pcoa <- function(data.in) {
  # Calculate convex hulls for each group
  hulls <- data.in %>%
    group_by(Cluster) %>%  # Use tidy evaluation for dynamic grouping
    slice(chull(NMDS1, NMDS2))           # Select points forming the convex hull
  
  # Create the PCA plot
  pcoa_plot <- ggplot(data.in) + 
    geom_point(aes(x = NMDS1, y = NMDS2, colour = Cluster), size = 4) + 
    geom_text(aes(x = NMDS1, y = NMDS2, label = Species), nudge_y = -0.01) +
    geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, 
                                   group = Cluster, 
                                   fill = Cluster), alpha = 0.3) +
    theme_bw() +
    xlab("NMDS1") + 
    ylab("NMDS2")
  
  return(pcoa_plot)
}

(reptile.pcoa = plot_pcoa(nmds_scores))
ggsave(plot = reptile.pcoa, filename = "plots/reptile.pcoa.plot.jpg", width = 6, height = 4, scale =2)

write.csv(reptile_traits, "data_clean/reptile_traits_clustered.csv")

#################
#### Frogs #### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
frog_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Amphibians") %>% mutate(stratum = as.factor(stratum),
                                                                                      nesting = as.factor(nesting),
                                                                                      diet= as.factor(diet),
                                                                                      dominant_pyrome = as.factor(dominant_pyrome))

# Traits of interest
all_traits <- c("scale(Mass_g)", "scale(biggest_patch_size)", "scale(n_habitat_patches)", "dominant_pyrome", "scale(pyrome_breadth)")
# Create the formula
response_var <- "meanben_firefreq"
smp_formula <- as.formula(paste(response_var, "~", paste(all_traits, collapse = " + ")))

# Fit the model as a Boosted Regression Tree
dat_modelling = frog_traits %>% drop_na(meanben_firefreq)

# Dataset is too small for BRTs

model = glm(traits_formula, data = dat_modelling)
all_traits <- c("Mass_g", "biggest_patch_size", "n_habitat_patches", "dominant_pyrome", "pyrome_breadth")

summary(model)

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
                                                                                       dominant_pyrome = as.factor(dominant_pyrome))
all_traits <- c("Mass_g", "biggest_patch_size", "n_habitat_patches", "dominant_pyrome", "pyrome_breadth")

traits_scaled = frog_traits %>% select(Common_Name,meanben_firefreq, all_of(all_traits)) %>% drop_na()

# Compute Gower distance matrix
gower_dist <- daisy(traits_scaled[,3:7], metric = "gower")

gmm_model <- Mclust(gower_dist, G=1:8) # Try between 1-20 clusters to see which is best. 

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
plot_pcoa <- function(data.in) {
  # Calculate convex hulls for each group
  hulls <- data.in %>%
    group_by(Cluster) %>%  # Use tidy evaluation for dynamic grouping
    slice(chull(NMDS1, NMDS2))           # Select points forming the convex hull
  
  # Create the PCA plot
  pcoa_plot <- ggplot(data.in) + 
    geom_point(aes(x = NMDS1, y = NMDS2, colour = Cluster), size = 4) + 
    geom_text(aes(x = NMDS1, y = NMDS2, label = Species), nudge_y = -0.01) +
    geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, 
                                   group = Cluster, 
                                   fill = Cluster), alpha = 0.3) +
    theme_bw() +
    xlab("NMDS1") + 
    ylab("NMDS2")
  
  return(pcoa_plot)
}

(frog.pcoa <- plot_pcoa(nmds_scores))

ggsave(plot = frog.pcoa, filename = "plots/frog.pcoa.plot.jpg", width = 6, height = 4, scale =2)

write.csv(frog_traits, "data_clean/frog_traits_clustered.csv")

#### ENDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

