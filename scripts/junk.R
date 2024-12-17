### Junk Script 


vic_fauna_traits_imputed = readRDS("data_clean/vic_fauna_traits_imputed.Rds")
dat = vic_fauna_traits_imputed

# Have a look at some metric comparisons
ggplot(dat, aes(x = meanben_firefreq, y = future_fire_impact)) + geom_point() + geom_smooth(method = "lm")
ggplot(dat, aes(x = meanben_firefreq, y = fame_lm_slope)) + geom_point() + geom_smooth(method = "lm")
ggplot(dat, aes(x = future_fire_impact, y = fame_lm_slope)) + geom_point() + geom_smooth(method = "lm")

# Define the trait groups
trait_groups = list(Movement_Physiology = c("scale(Mass_g)", "home_range_km2", "dispersal_km",  "max_longevity_d"),
                    Reproduction = c("litter_size_n", "litters_per_year_n","n_offspring_year"),
                    Habit = c("stratum", "stratum_aerial","stratum_arboreal_insessorial",
                              "stratum_aquatic", "stratum_fossorial", "stratum_generalist", "stratum_saxicolous",
                              "stratum_terrestrial", "volant"),
                    Behaviour = c("hibernation_torpor"),
                    Diet = c("diet_breadth_n","diet","diet_breadth_n", "diet_carnivore",
                             "diet_herbivore","diet_invertivore", "diet_infloresence",
                             "diet_omnivore","diet_granivore" ),
                    Biogeography = c("biggest_patch_size", "n_habitat_patches",
                                     "patch_isolation", "pyrome_breadth"))

mammals = dat %>% filter(Taxa_Group == "Mammals")
m = lm(meanben_firefreq ~ scale(Mass_g) + stratum + volant + nest_hollows + pyrome_breadth, data = mammals)
summary(m)

birds = dat %>% filter(Taxa_Group == "Birds")
m = lm(meanben_firefreq ~ scale(Mass_g) + stratum + diet + pyrome_breadth, data = birds)
plot(effects::allEffects(m))
summary(m)

reptiles = dat %>% filter(Taxa_Group == "Reptiles")
m = lm(meanben_firefreq ~ scale(Mass_g) +  diet + pyrome_breadth, data = reptiles)
plot(effects::allEffects(m))
summary(m)

frogs = dat %>% filter(Taxa_Group == "Amphibians")
m = lm(meanben_firefreq ~ scale(Mass_g) + pyrome_breadth, data = frogs)
plot(effects::allEffects(m))
summary(m)
?effects::allEffects

# Loop through each trait group
data.out = data.frame()
for (group in names(trait_groups)) {
  traits <- trait_groups[[group]]
  # Create a data frame to store AIC and adjusted R-squared for each model in the group
  group_results <- data.frame(Trait = character(length(traits)),
                              Group = character(length(traits)),
                              AIC = numeric(length(traits)),
                              Adj_R2 = numeric(length(traits)))
  # Fit models for each trait in the group
  for (i in seq_along(traits)) {
    trait <- traits[i]
    model <- lm(as.formula(paste("meanben_firefreq ~", trait)), data = dat)
    
    # Store the AIC and adjusted R-squared values
    group_results$Group[i] <- group
    group_results$Trait[i] <- trait
    group_results$AIC[i] <- AIC(model)
    group_results$Adj_R2_m[i] <- MuMIn::r.squaredGLMM(model)[1]
    group_results$Adj_R2_c[i] <- MuMIn::r.squaredGLMM(model)[2]
  }
  
  # Save the results for each group in the list
  data.out = rbind(data.out, group_results)
}

################# 
#### Mammals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

vic_mammal_traits_imputed = readRDS("data_clean/vic_mammal_traits_imputed.Rds")
dat = vic_mammal_traits_imputed

# Have a look at some metric comparisons
ggplot(dat, aes(x = meanBen_FireFreq, y = FutureImpact)) + geom_point() + geom_smooth(method = "lm")
ggplot(dat, aes(x = meanBen_FireFreq, y = lm_slope)) + geom_point() + geom_smooth(method = "lm")
ggplot(dat, aes(x = FutureImpact, y = lm_slope)) + geom_point() + geom_smooth(method = "lm")

dat$terrestrial_volant = as.numeric(dat$terrestrial_volant)
dat$foraging_stratum = as.factor(dat$foraging_stratum)
# Define the trait groups
trait_groups = list(Movement_Physiology = c("Mass_g", "home_range_km2", "dispersal_km",  "max_longevity_d"),
                    Reproduction = c("litter_size_n", "litters_per_year_n","n_offspring_year"),
                    Habitat = c("habitat_breadth_n"),
                    Habit = c("terrestrial_volant", "foraging_stratum", "fossoriality", "Arboreal_Insessorial", "Aerial", "Terrestrial"),
                    Behaviour = c("hibernation_torpor"),
                    Diet = c("det_diet_breadth_n","dphy_invertebrate","dphy_vertebrate", "dphy_plant"),
                    Biogeography = c("biggest_patch_size","n_patches","nned", "pyrome_breadth"))

# Loop through each trait group
data.out = data.frame()
for (group in names(trait_groups)) {
  traits <- trait_groups[[group]]
  # Create a data frame to store AIC and adjusted R-squared for each model in the group
  group_results <- data.frame(Trait = character(length(traits)),
                              Group = character(length(traits)),
                              AIC = numeric(length(traits)),
                              Adj_R2 = numeric(length(traits)))
  # Fit models for each trait in the group
  for (i in seq_along(traits)) {
    trait <- traits[i]
    model <- lm(as.formula(paste("meanBen_FireFreq ~", trait)), data = dat)
    
    # Store the AIC and adjusted R-squared values
    group_results$Group[i] <- group
    group_results$Trait[i] <- trait
    group_results$AIC[i] <- AIC(model)
    group_results$Adj_R2[i] <- summary(model)$adj.r.squared
  }
  
  # Save the results for each group in the list
  data.out = rbind(data.out, group_results)
}

best_trait_model = lm(meanBen_FireFreq ~ Mass_g + Arboreal_Insessorial + nned, data = dat)
summary(best_trait_model)

library(randomForest)

dat$Arboreal_Insessorial = as.factor(dat$Arboreal_Insessorial)
dat_rf = dat %>% drop_na(meanBen_FireFreq)
model <- randomForest(meanBen_FireFreq ~ Mass_g + dphy_plant + n_offspring_year + foraging_stratum + pyrome_breadth, data = dat_rf, importance = TRUE)

model
varImpPlot(model)  # Plot variable importance
summary(model)

partialPlot(model, pred.data = dat_rf, x.var = "Mass_g")
partialPlot(model, pred.data = dat_rf, x.var = "foraging_stratum")
partialPlot(model, pred.data = dat_rf, x.var = "pyrome_breadth")
partialPlot(model, pred.data = dat_rf, x.var = "dphy_plant")
partialPlot(model, pred.data = dat_rf, x.var = "n_offspring_year")

library(mclust)

dat_clust = dat_rf %>% 
  select(TAXON_ID, Sci_Name, 
         Mass_g, dphy_plant,foraging_stratum, n_offspring_year)

# Compute Gower distance matrix
library(cluster)

# Hierarchical clustering with Gower distance
gower_hc <- hclust(gower_dist, method = "ward.D2")

# Cut the dendrogram into clusters
clusters <- cutree(gower_hc, k = 4)
dat_clust$cluster <- clusters

# Plot dendrogram
plot(gower_hc)

library(vegan)
# Compute Gower distance matrix
gower_dist <- daisy(dat_clust[, 3:6], metric = "gower")
nmds <- metaMDS(gower_dist, k = 2)

# Extract NMDS scores and perform clustering
nmds_scores <- as.data.frame(scores(nmds))
nmds_kmeans <- kmeans(nmds_scores, centers = 6)
dat_clust$cluster <- nmds_kmeans$cluster

# Plot NMDS with clusters
plot(nmds_scores$NMDS1, nmds_scores$NMDS2, col = dat_clust$cluster, pch = 19)

# Perform PCA on standardized data
pca_model <- prcomp(dat_clust[,3:5], scale. = TRUE)

# Extract the first 2 or 3 principal components
pca_scores <- pca_model$x[, 1:3]

# Apply k-means on PCA scores
pca_kmeans <- kmeans(pca_scores, centers = 4, nstart = 25)
dat_clust$cluster <- pca_kmeans$cluster

# Plot PCA results with clusters
plot(pca_scores[, 1:2], col = dat_clust$cluster, pch = 19)

# Fit the model-based clustering
mclust_model <- Mclust(dat_clust[,3:6])

# Extract the optimal number of clusters and classification
dat_clust$cluster <- mclust_model$classification

# Visualize clusters
plot(mclust_model, what = "classification")

################# 
#### Birds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

vic_bird_traits_imputed = readRDS("data_clean/vic_bird_traits_imputed.Rds")
dat = vic_bird_traits_imputed

# Have a look at some metric comparisons
ggplot(dat, aes(x = meanBen_FireFreq, y = FutureImpact)) + geom_point() + geom_smooth(method = "lm")
ggplot(dat, aes(x = meanBen_FireFreq, y = lm_slope)) + geom_point() + geom_smooth(method = "lm")
ggplot(dat, aes(x = FutureImpact, y = lm_slope)) + geom_point() + geom_smooth(method = "lm")


# Define the trait groups
trait_groups = list(Movement_Physiology = c("Mass_g", "maximum_longevity_y"),
                    Reproduction = c("litters_or_clutches_per_y"),
                    Habitat = c("Habitat"),
                    Habit = c("Primary.Lifestyle", "Generalist", "Aquatic", "Arboreal_Insessorial", "Aerial", "Terrestrial"),
                    Diet = c("Trophic.Niche", "Omnivore","Vertivore","Invertivore","Granivore","Herbivore_terrestrial",
                             "Aquatic_predator","Nectarivore","Frugivore","Herbivore_aquatic"),
                    Biogeography = c("biggest_patch_size","n_patches","nned", "pyrome_breadth"))

# Loop through each trait group
data.out = data.frame()
for (group in names(trait_groups)) {
  traits <- trait_groups[[group]]
  # Create a data frame to store AIC and adjusted R-squared for each model in the group
  group_results <- data.frame(Trait = character(length(traits)),
                              Group = character(length(traits)),
                              AIC = numeric(length(traits)),
                              Adj_R2 = numeric(length(traits)))
  # Fit models for each trait in the group
  for (i in seq_along(traits)) {
    trait <- traits[i]
    model <- lm(as.formula(paste("meanBen_FireFreq ~", trait)), data = dat)
    
    # Store the AIC and adjusted R-squared values
    group_results$Group[i] <- group
    group_results$Trait[i] <- trait
    group_results$AIC[i] <- AIC(model)
    group_results$Adj_R2[i] <- summary(model)$adj.r.squared
  }
  
  # Save the results for each group in the list
  data.out = rbind(data.out, group_results)
}

best_trait_model = lm(meanBen_FireFreq ~ Mass_g + litters_or_clutches_per_y + 
                        Habitat + Trophic.Niche + Arboreal_Insessorial + pyrome_breadth, data = dat)
summary(best_trait_model)

dat_rf = dat %>% drop_na(meanBen_FireFreq, pyrome_breadth) %>%
  mutate(Habitat = as.factor(Habitat), Trophic.Niche = as.factor(Trophic.Niche))

model <- randomForest(meanBen_FireFreq ~ Mass_g + litters_or_clutches_per_y + 
                        Habitat + Trophic.Niche + Arboreal_Insessorial + pyrome_breadth, data = dat_rf, importance = TRUE)

varImpPlot(model)  # Plot variable importance
summary(model)

partialPlot(model, pred.data = dat_rf, x.var = "Mass_g")
partialPlot(model, pred.data = dat_rf, x.var = "litters_or_clutches_per_y")
partialPlot(model, pred.data = dat_rf, x.var = "Habitat")
partialPlot(model, pred.data = dat_rf, x.var = "Trophic.Niche")
partialPlot(model, pred.data = dat_rf, x.var = "Arboreal_Insessorial")
partialPlot(model, pred.data = dat_rf, x.var = "pyrome_breadth")


################# 
#### Reptiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
vic_reptile_traits_imputed = readRDS("data_clean/vic_reptile_traits_imputed.Rds")
dat = vic_reptile_traits_imputed
dat$litter_size_n = as.numeric(dat$litter_size_n)
# Define the trait groups
trait_groups = list(Movement_Physiology = c("Mass_g"),
                    Reproduction = c("maximum_longevity_y","litter_size_n","litters_per_year_n"),
                    Habitat = c("Habitat"),
                    Habit = c("Arboreal_Insessorial","Saxicolous","Fossorial","Aquatic","Cryptic", "Terrestrial"),
                    Diet = c("Omnivore", "Carnivore", "Herbivore"),
                    Biogeography = c("biggest_patch_size","n_patches","nned", "pyrome_breadth"))

# Loop through each trait group
data.out = data.frame()
for (group in names(trait_groups)) {
  traits <- trait_groups[[group]]
  # Create a data frame to store AIC and adjusted R-squared for each model in the group
  group_results <- data.frame(Trait = character(length(traits)),
                              Group = character(length(traits)),
                              AIC = numeric(length(traits)),
                              Adj_R2 = numeric(length(traits)))
  # Fit models for each trait in the group
  for (i in seq_along(traits)) {
    trait <- traits[i]
    model <- lm(as.formula(paste("meanBen_FireFreq ~", trait)), data = dat)
    
    # Store the AIC and adjusted R-squared values
    group_results$Group[i] <- group
    group_results$Trait[i] <- trait
    group_results$AIC[i] <- AIC(model)
    group_results$Adj_R2[i] <- summary(model)$adj.r.squared
  }
  
  # Save the results for each group in the list
  data.out = rbind(data.out, group_results)
}

best_trait_model = lm(meanBen_FireFreq ~ litter_size_n + pyrome_breadth + Herbivore + Fossorial + Mass_g, data=dat)
summary(best_trait_model)

################# 
#### Frogs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
vic_frog_traits_imputed = readRDS("data_clean/vic_frog_traits_imputed.Rds")
dat = vic_frog_traits_imputed

# Define the trait groups
trait_groups = list(Movement_Physiology = c("Mass_g"),
                    Reproduction = c("Age_at_maturity_min_y", "maximum_longevity_y", "Reproductive_output_y"),
                    #Habitat = c("Habitat"),
                    Habit = c("Fossorial","Terrestrial","Aquatic", "Arboreal_Insessorial"),
                    #Diet = c("Omnivore", "Carnivore", "Herbivore"),
                    Biogeography = c("biggest_patch_size","n_patches","nned", "pyrome_breadth"))
# Loop through each trait group
data.out = data.frame()
for (group in names(trait_groups)) {
  traits <- trait_groups[[group]]
  # Create a data frame to store AIC and adjusted R-squared for each model in the group
  group_results <- data.frame(Trait = character(length(traits)),
                              Group = character(length(traits)),
                              AIC = numeric(length(traits)),
                              Adj_R2 = numeric(length(traits)))
  # Fit models for each trait in the group
  for (i in seq_along(traits)) {
    trait <- traits[i]
    model <- lm(as.formula(paste("meanBen_FireFreq ~", trait)), data = dat)
    
    # Store the AIC and adjusted R-squared values
    group_results$Group[i] <- group
    group_results$Trait[i] <- trait
    group_results$AIC[i] <- AIC(model)
    group_results$Adj_R2[i] <- summary(model)$adj.r.squared
  }
  
  # Save the results for each group in the list
  data.out = rbind(data.out, group_results)
}

best_trait_model = lm(meanBen_FireFreq ~ biggest_patch_size + Age_at_maturity_min_y + Arboreal_Insessorial + Mass_g, data=dat)
summary(best_trait_model)


vic_fauna_traits= readRDS("data_clean/vic_fauna_traits_imputed.Rds")
dat = vic_fauna_traits

# Have a look at some metric comparisons
ggplot(dat, aes(x = meanben_firefreq, y = future_fire_impact)) + geom_point() + geom_smooth(method = "lm")
ggplot(dat, aes(x = meanben_firefreq, y = fame_lm_slope)) + geom_point() + geom_smooth(method = "lm")
ggplot(dat, aes(x = future_fire_impact, y = fame_lm_slope)) + geom_point() + geom_smooth(method = "lm")

# Define the trait groups
trait_groups = list(Movement_Physiology = c("Mass_g", "home_range_km2", "dispersal_km",  "max_longevity_d"),
                    Reproduction = c("litter_size_n", "litters_per_year_n","n_offspring_year"),
                    Habit = c("stratum", "stratum_aerial","stratum_arboreal_insessorial",
                              "stratum_aquatic", "stratum_fossorial", "stratum_generalist", "stratum_saxicolous",
                              "stratum_terrestrial", "volant"),
                    Nesting = c("nesting","nest_burrow","nest_cave","nest_ground","nest_hollows",                
                                "nest_branch"),
                    Behaviour = c("hibernation_torpor"),
                    Diet = c("diet_breadth_n","diet","diet_breadth_n", "diet_carnivore",
                             "diet_herbivore","diet_invertivore", "diet_infloresence",
                             "diet_omnivore","diet_granivore"),
                    Biogeography = c("biggest_patch_size", "n_habitat_patches",
                                     "patch_isolation", "pyrome_breadth"))


# Loop through each trait group
data.out = data.frame()
for (group in names(trait_groups)) {
  traits <- trait_groups[[group]]
  # Create a data frame to store AIC and adjusted R-squared for each model in the group
  group_results <- data.frame(Trait = character(length(traits)),
                              Group = character(length(traits)),
                              sample_size = numeric(length(traits)),
                              AIC = numeric(length(traits)),
                              Adj_R2 = numeric(length(traits)))
  # Fit models for each trait in the group
  for (i in seq_along(traits)) {
    trait <- traits[i]
    model <- lm(as.formula(paste("meanben_firefreq ~", trait)), data = dat)
    
    # Store the AIC and adjusted R-squared values
    group_results$Group[i] <- group
    group_results$Trait[i] <- trait
    group_results$sample_size <- model$df.residual
    group_results$AIC[i] <- AIC(model)
    group_results$Adj_R2_m[i] <- MuMIn::r.squaredGLMM(model)[1]
    group_results$Adj_R2_c[i] <- MuMIn::r.squaredGLMM(model)[2]
  }
  
  # Save the results for each group in the list
  data.out = rbind(data.out, group_results)
}


dat_rf = dat %>% 
  drop_na(meanben_firefreq, Mass_g, nest_hollows, litter_size_n, pyrome_breadth) %>% 
  mutate(Taxa_Group = as.factor(Taxa_Group))

library(randomForest)

model <- randomForest(meanben_firefreq ~ 
                        Taxa_Group + stratum_aerial + stratum_arboreal_insessorial + Mass_g + pyrome_breadth, 
                      data = dat_rf, importance=TRUE)

summary(model)
model
varImpPlot(model)  # Plot variable importance
summary(model)

partialPlot(model, pred.data = dat_rf, x.var = "Mass_g")
partialPlot(model, pred.data = dat_rf, x.var = "Taxa_Group")
partialPlot(model, pred.data = dat_rf, x.var = "stratum_aerial")
partialPlot(model, pred.data = dat_rf, x.var = "pyrome_breadth")
partialPlot(model, pred.data = dat_rf, x.var = "litter_size_n")
partialPlot(model, pred.data = dat_rf, x.var = "n_offspring_year")

# Compute Gower distance matrix
library(cluster)
library(vegan)
dat_clust=dat_rf 
dat_clust = dat_clust %>% select(Mass_g, stratum_aerial, stratum_arboreal_insessorial, stratum_aquatic, stratum_cryptic, stratum_fossorial, stratum_generalist, stratum_saxicolous, 
                                 stratum_terrestrial, nest_burrow, nest_cave, nest_ground, nest_hollows, nest_branch, volant, hibernation_torpor, diet_carnivore, diet_herbivore, 
                                 diet_invertivore, diet_infloresence, diet_omnivore, diet_granivore, litter_size_n, biggest_patch_size, n_habitat_patches, patch_isolation, 
                                 pyrome_breadth, meanben_firefreq)
gower_dist <- daisy(dat_clust, metric = "gower")
nmds <- metaMDS(gower_dist,method="gower", k = 2, trymax = 30)
plot(nmds)
# Extract NMDS scores and perform clustering
nmds_scores <- as.data.frame(scores(nmds))

wss <- sapply(1:30, function(k) {
  kmeans(nmds_scores, centers = k, nstart = 25)$tot.withinss
})

plot(1:30, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares")


silhouette_scores <- sapply(2:30, function(k) {
  km <- kmeans(nmds_scores, centers = k, nstart = 25)
  ss <- silhouette(km$cluster, dist(nmds_scores))
  mean(ss[, 3])  # Average silhouette width
})

plot(2:30, silhouette_scores, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Average silhouette width")

# Gap statistic with 10 bootstraps
set.seed(123)
gap_stat <- clusGap(nmds_scores, FUN = kmeans, nstart = 25, K.max = 30, B = 10)
print(gap_stat, method = "firstmax")

# Plot the gap statistic
plot(gap_stat, frame = FALSE)

# Model based approach to determining the number of clusters from the data... 
library(NbClust)
set.seed(123)
nb <- NbClust(nmds_scores, distance = "euclidean", min.nc = 2, max.nc = 30, method = "kmeans")


nmds_kmeans <- kmeans(nmds_scores, centers = 20)
dat_clust=dat_rf
dat_clust$cluster <- nmds_kmeans$cluster
clusters = dat_clust %>% select(Taxa_Group, Common_Name, cluster)

# Plot NMDS with clusters
plot(nmds_scores$NMDS1, nmds_scores$NMDS2, col = dat_clust$cluster, pch = 19)

library(mboost)

glm1 <- glmboost(meanben_firefreq ~ scale(pyrome_breadth) + scale(Mass_g) + scale(nest_hollows), data = dat)
coef(glm1, off2int=TRUE)
plot(glm1, off2int = TRUE)
summary(glm1)


library(mclust)
# Fit the model-based clustering
dat_clust = dat_rf %>% 
  select(Taxon_ID, Scientific_Name, Taxa_Group, 
         Mass_g, Taxa_Group, stratum_aerial, stratum_arboreal_insessorial, Mass_g, pyrome_breadth)
mclust_model <- Mclust(dat_clust)

# Extract the optimal number of clusters and classification
dat_clust$cluster <- mclust_model$classification





# 
# #### Mammals
# # Amniote and otherone basically the same. Use Amniote
# mammal.traits = read.csv("data_raw/Amniote_Database_Aug_2015.csv")
# 
# mammal.traits = mammal.traits %>% mutate(Binomial = paste(genus, species)) %>%
#   select(Binomial, litter_or_clutch_size_n) %>% 
#   rename("Litter_Size" = "litter_or_clutch_size_n") %>% 
#   mutate(Litter_Size = ifelse(Litter_Size == -999, NA, Litter_Size))
# 
# mammals = left_join(mammals, mammal.traits, by=c("SCIENTIFIC_NAME" = "Binomial"))
# 
# #### Birds
# # Amniote is better
# bird.traits = read.csv("data_raw/Amniote_Database_Aug_2015.csv")
# 
# bird.traits = bird.traits %>% mutate(Binomial = paste(genus, species)) %>%
#   select(Binomial, litter_or_clutch_size_n) %>%
#   mutate(Binomial = gsub("Puffinus", "Ardenna", Binomial)) %>%
#   rename("Litter_Size" = "litter_or_clutch_size_n") %>% 
#   mutate(Litter_Size = ifelse(Litter_Size == -999, NA, Litter_Size))
# 
# birds = left_join(birds, bird.traits, by=c("SCIENTIFIC_NAME" = "Binomial"))
# 
# missing.birds = df_birds %>% filter(is.na(Litter_Size)==TRUE) %>% select(-Litter_Size)
# missing.birds.data = left_join(missing.birds, bird.traits, by=c("species" = "Binomial")) %>% group_by(species_label) %>% summarise(Litter_Size = mean(Litter_Size,na.rm=TRUE))
# missing.birds = left_join(missing.birds, missing.birds.data, by="species_label")
# missing.birds$Litter_Size <- ifelse(missing.birds$species_label=="Petrochelidon ariel", 4, missing.birds$Litter_Size) # Species from same genus
# missing.birds$Litter_Size <- ifelse(missing.birds$species_label=="Petrochelidon nigricans", 4, missing.birds$Litter_Size) # Species from same genus
# missing.birds$Litter_Size <- ifelse(missing.birds$species_label=="Taeniopygia castanotis", 4.5, missing.birds$Litter_Size) # Species from same genus (Zebra Finch)
# 
# df_birds = df_birds %>% drop_na(Litter_Size) %>% rbind(missing.birds)
# 
# #### Reptiles
# reptile.traits = read.csv("data_raw/Amniote_Database_Aug_2015.csv")
# 
# reptile.traits = reptile.traits %>% mutate(Binomial = paste(genus, species)) %>%
#   select(Binomial, litter_or_clutch_size_n) %>% 
#   rename("Litter_Size" = "litter_or_clutch_size_n")  %>% 
#   mutate(Litter_Size = ifelse(Litter_Size == -999, NA, Litter_Size))
# 
# reptiles = left_join(reptiles, reptile.traits, by=c("SCIENTIFIC_NAME" = "Binomial"))
# 
# #### Frogs
# frog.traits = read.csv("data_raw/amphibians/AmphiBIO_v1.csv") 
# frog.traits = frog.traits %>% 
#   select(Species,Litter_size_max_n) %>% rename("Litter_Size" = "Litter_size_max_n")
# 
# frogs = left_join(frogs, frog.traits, by = c("SCIENTIFIC_NAME"="Species"))
# 
# #### Diet ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #### Mammals
# mammal.diet = read.table("data_raw/traits/mammals/MamFuncDat.txt", header=TRUE, sep= "\t")
# mammal.traits = mammal.diet %>% select(Scientific, 
#                                        Diet.Inv,    # Invertebrate
#                                        Diet.Vend,   # birds and Birds
#                                        Diet.Vect,   # Reptiles and amphibians
#                                        Diet.Vfish,   # Fish
#                                        Diet.Vunk,   # Vetebrates General or Unknown
#                                        Diet.Scav,   # Scavenger
#                                        Diet.Fruit,  # Fruit
#                                        Diet.Nect,   # Nectar and pollen
#                                        Diet.Seed,   # Seed
#                                        Diet.PlantO) %>% # Other Plant Matter
#   mutate(D_Vertebrates = rowSums(.[, c("Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk", "Diet.Scav")], na.rm = TRUE),
#          D_Invertebrates = Diet.Inv, 
#          D_Fruits = Diet.Fruit, 
#          D_Seeds = Diet.Seed, 
#          D_Flowers = Diet.Nect,
#          D_PlantO = Diet.PlantO) %>% select(Scientific, D_Vertebrates, D_Invertebrates, D_Fruits, D_Seeds, D_Flowers, D_PlantO)
# 
# df_mammals = left_join(df_mammals, mammal.traits, by =c("species_label"="Scientific")) 
# 
# 
# #### Birds
# bird.diet = read.csv("data_raw/traits/birds/BirdFuncDat.csv")
# 
# bird.traits = bird.diet %>% select(Scientific, 
#                                    Diet.Inv,    # Invertebrate
#                                    Diet.Vend,   # birds and Birds
#                                    Diet.Vect,   # Reptiles and amphibians
#                                    Diet.Vfish,   # Fish
#                                    Diet.Vunk,   # Vetebrates General or Unknown
#                                    Diet.Scav,   # Scavenger
#                                    Diet.Fruit,  # Fruit
#                                    Diet.Nect,   # Nectar and pollen
#                                    Diet.Seed,   # Seed
#                                    Diet.PlantO) %>% # Other Plant Matter
#   mutate(D_Vertebrates = rowSums(.[, c("Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk", "Diet.Scav")], na.rm = TRUE),
#          D_Invertebrates = Diet.Inv, 
#          D_Fruits = Diet.Fruit, 
#          D_Seeds = Diet.Seed, 
#          D_Flowers = Diet.Nect,
#          D_PlantO = Diet.PlantO) %>% select(Scientific, D_Vertebrates, D_Invertebrates, D_Fruits, D_Seeds, D_Flowers, D_PlantO)
# 
# df_birds = left_join(df_birds, bird.traits, by =c("species_label"="Scientific"))
# 
# # Missing species
# missing.species = df_birds %>% filter(is.na(df_birds$D_Vertebrates)==TRUE) %>% 
#   select(-c(D_Vertebrates, D_Invertebrates, D_Fruits, D_Seeds, D_Flowers, D_PlantO))
# 
# missing.birds.data = left_join(missing.species, bird.traits, by=c("species" = "Scientific")) %>% 
#   group_by(species_label) %>% summarise(D_Vertebrates = mean(D_Vertebrates, na.rm=TRUE),
#                                         D_Invertebrates = mean(D_Invertebrates, na.rm=TRUE),
#                                         D_Fruits = mean(D_Fruits, na.rm=TRUE),
#                                         D_Seeds = mean(D_Seeds, na.rm=TRUE),
#                                         D_Flowers = mean(D_Flowers, na.rm = TRUE),
#                                         D_PlantO = mean(D_PlantO, na.rm=TRUE))
# # Missing species part two 
# missing.birds2 = missing.birds.data %>% 
#   filter(is.na(missing.birds.data$D_Vertebrates)) %>% 
#   select(-c(D_Vertebrates, D_Invertebrates, D_Fruits, D_Seeds, D_Flowers, D_PlantO))
# 
# # These are basically all just changes to taxonomy.. but same species. A couple needed to use close relatives. 
# missing.birds2$Ref_Species = c("Xanthomyza phrygia", 
#                                "Chalcophaps indica", 
#                                "Cacatua roseicapilla",
#                                "Phylidonyris melanops",
#                                "Cuculus pallidus",
#                                "Hirundo ariel", 
#                                "Hirundo nigricans",
#                                "Petroica multicolor",
#                                "Chthonicola sagittatus",
#                                "Rhipidura fuliginosa", 
#                                "Stigmatopelia chinensis",
#                                "Monarcha trivirgatus", 
#                                "Taeniopygia guttata",
#                                "Trichoglossus haematodus")
# 
# missing.birds2 = left_join(missing.birds2, bird.traits, by=c("Ref_Species" = "Scientific")) %>% 
#   group_by(species_label) %>% summarise(D_Vertebrates = mean(D_Vertebrates, na.rm=TRUE),
#                                         D_Invertebrates = mean(D_Invertebrates, na.rm=TRUE),
#                                         D_Fruits = mean(D_Fruits, na.rm=TRUE),
#                                         D_Seeds = mean(D_Seeds, na.rm=TRUE),
#                                         D_Flowers = mean(D_Flowers, na.rm = TRUE),
#                                         D_PlantO = mean(D_PlantO, na.rm=TRUE))
# 
# missing.birds.data = drop_na(missing.birds.data, D_Vertebrates)
# missing.birds.data = rbind(missing.birds.data, missing.birds2)
# 
# missing.birds = left_join(missing.species, missing.birds.data, by="species_label")
# 
# df_birds = df_birds %>% drop_na(D_Vertebrates) %>% rbind(missing.birds)
# 
# 
# #### Reptiles
# reptile.traits = readxl::read_excel("data_raw/traits/reptiles/ReptTraits%20dataset%20v1-1.xlsx") %>%
#   select(Species, Diet, `Diet: comments`)
# 
# df_reptiles = filter(species_df, TAX_GROUP == "reptiles")
# df_reptiles = left_join(df_reptiles, reptile.traits, by =c("spp_label" = "Species"))
# 
# #### Frogs
# frog.traits = read.csv("data_raw/traits/amphibians/AmphiBIO_v1.csv") 
# frog.traits = frog.traits %>% 
#   select(Species, 
#          Vert, 
#          Arthro,
#          Leaves, 
#          Flowers, 
#          Seeds, 
#          Fruits) %>% rename("D_Invertebrates" = "Arthro", "D_Vertebrates" = "Vert", "D_Fruits" = "Fruits", 
#                             "D_Seeds" = "Seeds", "D_Flowers" = "Flowers", "D_PlantO" = "Leaves") %>%
#   select(Species, D_Vertebrates, D_Invertebrates, D_Fruits, D_Seeds, D_Flowers, D_PlantO)
# 
# df_frogs = left_join(df_frogs, frog.traits, by = c("species_label"="Species"))
# 
# #### Microhabitat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #### Mammals
# m.traits = read.csv("data_raw/traits/mammals/combine_trait_data_imputed.csv")
# mammal.traits = m.traits %>% select(iucn2020_binomial, fossoriality, foraging_stratum) %>%
#   mutate(Fossorial = ifelse(fossoriality == 1, 1,
#                             ifelse(fossoriality ==2,0, NA))) %>%
#   mutate(Aquatic = ifelse(foraging_stratum == "M", "1",
#                           ifelse(is.na(foraging_stratum)==TRUE, NA, 0)),
#          Terrestrial = ifelse(foraging_stratum == "G", "1",
#                               ifelse(is.na(foraging_stratum)==TRUE, NA, 0)),
#          Arboreal = ifelse(foraging_stratum %in% c("Ar", "S"), "1",
#                            ifelse(is.na(foraging_stratum)==TRUE, NA, 0)),
#          Aerial = ifelse(foraging_stratum == "A", "1",
#                          ifelse(is.na(foraging_stratum)==TRUE, NA, 0))) %>%
#   select(-c(fossoriality, foraging_stratum))
# 
# df_mammals = left_join(df_mammals, mammal.traits, by=c("species_label" = "iucn2020_binomial"))
# 
# missing.mammals = df_mammals %>% 
#   filter(is.na(df_mammals$Fossorial)) %>% 
#   select(X:D_PlantO)
# 
# missing.mammals$Ref_species = c("Antechinus agilis", "Canis latrans", "Felis silvestris", "Austronomus australis")
# missing.mammals = left_join(missing.mammals, mammal.traits, by=c("Ref_species" = "iucn2020_binomial")) %>% select(-Ref_species)
# df_mammals = df_mammals %>% drop_na(Fossorial) %>% rbind(missing.mammals)
# 
# #### Birds
# birds.traits = read.csv("data_raw/traits/birds/AVONET1_BirdLife.csv")
# birds.traits = birds.traits %>% 
#   select(Species1,Primary.Lifestyle) %>%
#   mutate(Aquatic = ifelse(Primary.Lifestyle == "Aquatic", 1,0),
#          Terrestrial = ifelse(Primary.Lifestyle == "Terrestrial", 1,0),
#          Arboreal = ifelse(Primary.Lifestyle == "Insessorial", 1,0),
#          Aerial = ifelse(Primary.Lifestyle == "Aerial", 1,0),
#          Fossorial = 0) 
# 
# df_birds = left_join(df_birds, birds.traits, by = c("species_label" = "Species1")) %>%
#   mutate(Terrestrial = ifelse(Primary.Lifestyle=="Generalist" & Aquatic == 0, 1, Terrestrial),
#          Arboreal = ifelse(Primary.Lifestyle=="Generalist" & Aquatic == 0, 1, Arboreal),
#          Aerial = ifelse(Primary.Lifestyle=="Generalist" & Aquatic == 0, 1, Aerial),) %>%
#   select(-Primary.Lifestyle)
# 
# #### Reptiles
# reptile.traits =  read_excel("data_raw/traits/reptiles/ReptTraits%20dataset%20v1-1.xlsx")
# reptile.traits = reptile.traits %>% 
#   select(Species,Microhabitat) %>%
#   mutate(Aquatic = ifelse(str_detect(Microhabitat, "Aquatic|Semiaquatic|Marine"), 1, ifelse(Microhabitat=="NA", NA, 0)),
#          Terrestrial = ifelse(str_detect(Microhabitat, "Terrestrial|Cryptic|Saxicolous"), 1, ifelse(Microhabitat=="NA", NA, 0)),
#          Arboreal = ifelse(str_detect(Microhabitat, "Arboreal"), 1, ifelse(Microhabitat=="NA", NA, 0)),
#          Aerial = ifelse(Microhabitat == "Aerial", 1, ifelse(Microhabitat=="NA", NA, 0)),
#          Fossorial = ifelse(str_detect(Microhabitat, "Fossorial"), 1, ifelse(Microhabitat=="NA", NA, 0))) %>%
#   select(-Microhabitat)
# 
# df_reptiles = left_join(df_reptiles, reptile.traits, by = c("species_label"="Species"))
# 
# #### Frogs
# frog.traits = read.csv("data_raw/traits/amphibians/AmphiBIO_v1.csv") 
# 
# frog.traits = frog.traits %>% 
#   select(Species, Fos, Ter, Aqu, Arb) %>%
#   mutate(Aquatic = ifelse(Aqu==1, 1, 0),
#          Fossorial = ifelse(Fos == 1, 1, 0),
#          Terrestrial = ifelse(Ter==1, 1, 0), 
#          Arboreal = ifelse(Arb ==1, 1, 0),
#          Aerial = 0) %>% select(-Fos, -Ter, -Aqu, -Arb) %>%
#   mutate_at(vars(Aquatic:Aerial), ~ifelse(is.na(.), 0, .))
# 
# df_frogs = left_join(df_frogs, frog.traits, by = c("species_label"="Species"))