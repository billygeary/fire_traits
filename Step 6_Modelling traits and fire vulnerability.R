## Trait modelling

#################
#### Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
library(tidyverse)
library(brms)
library(ggplot2)

# Read in the imputed data
vic_fauna_traits_imputed = readRDS("data_clean/vic_fauna_traits_imputed.Rds")

#################
#### Mammals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
mammal_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Mammals")

# Define the trait groups
trait_groups = list(Movement_Physiology = c("Mass_g", "home_range_km2", "dispersal_km",  "max_longevity_d"),
                    Reproduction = c("litter_size_n", "litters_per_year_n","n_offspring_year"),
                    Habit = c("stratum","stratum_aerial","stratum_arboreal_insessorial", "stratum_aquatic","stratum_cryptic","stratum_fossorial"),
                    Nesting = c("nesting", "nest_burrow","nest_cave","nest_ground","nest_hollows","nest_branch"),
                    Behaviour = c("hibernation_torpor", "volant"),
                    Diet = c("diet","diet_breadth_n","diet_animals","diet_plants", "diet_carnivore","diet_herbivore",
                             "diet_invertivore","diet_infloresence","diet_omnivore","diet_granivore"),
                    Biogeography = c("biggest_patch_size","n_habitat_patches","patch_isolation", "pyrome_breadth"))

dat = mammal_traits

# Loop through each trait group
data.out <- data.frame()

# Loop through each trait group
for (group in names(trait_groups)) {
  traits <- trait_groups[[group]]
  
  # Create a data frame to store AIC, adjusted R-squared, and errors for each model in the group
  group_results <- data.frame(Trait = character(length(traits)),
                              Group = character(length(traits)),
                              AIC = numeric(length(traits)),
                              Adj_R2 = numeric(length(traits)),
                              Error = character(length(traits)),
                              stringsAsFactors = FALSE)
  
  # Fit models for each trait in the group
  for (i in seq_along(traits)) {
    trait <- traits[i]
    
    # Use tryCatch to handle errors
    result <- tryCatch({
      model <- lm(as.formula(paste("meanben_firefreq ~", trait)), data = dat)
      
      # Store AIC and adjusted R-squared values
      list(AIC = AIC(model), Adj_R2 = summary(model)$adj.r.squared, Error = NA)
    }, error = function(e) {
      # Capture error message
      list(AIC = NA, Adj_R2 = NA, Error = as.character(e$message))
    })
    
    # Store results in the group_results data frame
    group_results$Group[i] <- group
    group_results$Trait[i] <- trait
    group_results$AIC[i] <- result$AIC
    group_results$Adj_R2[i] <- result$Adj_R2
    group_results$Error[i] <- result$Error
  }
  
  # Append the group's results to the main data.out data frame
  data.out <- rbind(data.out, group_results)
}

# View final results
data.out

# Fit a model
dat = dat %>% drop_na(meanben_firefreq)
options(na.action = "na.omit")
mammals = lm(meanben_firefreq ~ scale(Mass_g) + scale(litter_size_n) + stratum + nesting + hibernation_torpor + volant + 
               diet_animals + diet_plants + scale(patch_isolation) + scale(pyrome_breadth), data = dat)

options(na.action = "na.fail")
mammals.dredge = dredge(mammals)

library(brms)
mammal.best = brm(meanben_firefreq ~ scale(Mass_g) + scale(patch_isolation) + scale(pyrome_breadth) + 
                    stratum + nesting, data = dat)

summary(mammal.best)

conditional_effects(mammal.best, "nesting")
conditional_effects(mammal.best, "stratum")
conditional_effects(mammal.best, "patch_isolation")
conditional_effects(mammal.best, "pyrome_breadth")
conditional_effects(mammal.best, "Mass_g")

#################
#### Birds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
birds_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Birds")
options(na.action = "na.omit")
dat = birds_traits

# Loop through each trait group
data.out <- data.frame()

# Loop through each trait group
for (group in names(trait_groups)) {
  traits <- trait_groups[[group]]
  
  # Create a data frame to store AIC, adjusted R-squared, and errors for each model in the group
  group_results <- data.frame(Trait = character(length(traits)),
                              Group = character(length(traits)),
                              AIC = numeric(length(traits)),
                              Adj_R2 = numeric(length(traits)),
                              Error = character(length(traits)),
                              stringsAsFactors = FALSE)
  
  # Fit models for each trait in the group
  for (i in seq_along(traits)) {
    trait <- traits[i]
    
    # Use tryCatch to handle errors
    result <- tryCatch({
      model <- lm(as.formula(paste("meanben_firefreq ~", trait)), data = dat)
      
      # Store AIC and adjusted R-squared values
      list(AIC = AIC(model), Adj_R2 = summary(model)$adj.r.squared, Error = NA)
    }, error = function(e) {
      # Capture error message
      list(AIC = NA, Adj_R2 = NA, Error = as.character(e$message))
    })
    
    # Store results in the group_results data frame
    group_results$Group[i] <- group
    group_results$Trait[i] <- trait
    group_results$AIC[i] <- result$AIC
    group_results$Adj_R2[i] <- result$Adj_R2
    group_results$Error[i] <- result$Error
  }
  
  # Append the group's results to the main data.out data frame
  data.out <- rbind(data.out, group_results)
}

# Fit a model
dat = dat %>% drop_na(meanben_firefreq)
options(na.action = "na.omit")
birds = brm(meanben_firefreq ~ scale(Mass_g) + scale(litters_per_year_n) + stratum + nesting + diet + 
               scale(patch_isolation) + scale(biggest_patch_size) + scale(pyrome_breadth), data = dat)

brms::conditional_effects(birds, "biggest_patch_size")

options(na.action = "na.fail")
birds.dredge = dredge(birds)

head(birds.dredge)

birds.best = lm(meanben_firefreq ~ scale(Mass_g) + scale(patch_isolation) + scale(pyrome_breadth) + stratum, dat = dat)
summary(birds.best)
plot(effects::allEffects(birds.best))

#################
#### Reptiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
mammal_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Mammals")


#################
#### Frogs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
mammal_traits = vic_fauna_traits_imputed %>% filter(Taxa_Group=="Mammals")
