# ------------------------------------------------------------------------------
# Script Name: Step 5_Trait Imputation
# Author: Billy Geary
# Date Created: 17/12/2024
# Last Modified: 17/12/2024
# Purpose: This script uses Trait Imputation methods to fill in the gaps in the trait information gathered in Step 4. Trait imputation is conducted for each taxon
#          group separately, using the MICE package. 
# Outputs: A database of the same structure as teh one created in Step 4, but with the imputed trait information also included.   
# ------------------------------------------------------------------------------

## Set up
library(mice)
library(tidyverse)
#################
#### All Fauna ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
vic_fauna_traits= readRDS("data_clean/vic_fauna_traits.Rds")

vic_fauna_traits = vic_fauna_traits %>% mutate(stratum = as.factor(stratum),
                                            nesting = as.factor(nesting),
                                            diet = as.factor(diet),
                                            dominant_pyrome = as.factor(dominant_pyrome))

pred.mat = as.matrix(c(0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
pred.mat = matrix(rep(pred.mat, each = length(pred.mat)), nrow = length(pred.mat), byrow = TRUE)

################# 
#### Mammals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
vic_mammal_traits = vic_fauna_traits %>% filter(Taxa_Group=="Mammals") %>% droplevels()

# Define the column(s) you do not want to impute
# Create a default method vector using `mice` defaults and set `""` for excluded columns
methods <- make.method(vic_mammal_traits[,1:40], defaultMethod = c("rf", "rf", "polyreg", "polr"))

imputed_data <- mice(vic_mammal_traits[,1:40], 
                     m = 5, # How many imputations
                     predictorMatrix = pred.mat,
                     method = methods, # Method used (needs to be specific to continuous vs. categorical)
                     maxit = 50, # Number of iterations for algorithm 
                     seed = 123) # Set Seed

# Bring imputed data back into the dataset 
vic_mammal_traits_imputed <- complete(imputed_data)
vic_mammal_traits_imputed <- cbind(vic_mammal_traits_imputed, vic_mammal_traits[,41:49])

# Save
saveRDS(vic_mammal_traits_imputed, "data_clean/vic_mammal_traits_imputed.Rds")# Mammals

################# 
#### Birds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
vic_bird_traits = vic_fauna_traits %>% filter(Taxa_Group=="Birds") %>% droplevels()

#vic_bird_traits = readRDS("data_clean/vic_birds_traits.Rds")

mice::md.pattern(vic_bird_traits)

# Do the imputation
names(vic_bird_traits)
#pred.mat = as.matrix(c(0,0,0,0,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
#pred.mat = matrix(rep(pred.mat, each = length(pred.mat)), nrow = length(pred.mat), byrow = TRUE)

# Define the column(s) you do not want to impute
# Create a default method vector using `mice` defaults and set `""` for excluded columns
methods <- make.method(vic_bird_traits[,1:40], defaultMethod = c("rf", "rf", "polyreg", "polr"))

imputed_data <- mice(vic_bird_traits[,1:40], 
                     m = 5, # How many imputations
                     predictorMatrix = pred.mat,
                     method = methods, # Method used (needs to be specific to continuous vs. categorical)
                     maxit = 50, # Number of iterations for algorithm 
                     seed = 123) # Set Seed

# Bring imputed data back into the dataset 
vic_bird_traits_imputed <- complete(imputed_data)
vic_bird_traits_imputed <- cbind(vic_bird_traits_imputed, vic_bird_traits[,41:49])

# Save
saveRDS(vic_bird_traits_imputed, "data_clean/vic_bird_traits_imputed.Rds")

################# 
#### Reptiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
vic_reptile_traits = vic_fauna_traits %>% filter(Taxa_Group=="Reptiles") %>% droplevels()

#vic_reptile_traits = readRDS("data_clean/vic_reptile_traits.Rds")

mice::md.pattern(vic_reptile_traits)

# Do the imputation
names(vic_reptile_traits)

# Pick the traits we actually have enough data to impute - the rest will stay NAs
vic_reptile_traits_to_impute = vic_reptile_traits %>% select(c(1:11,13:20,29,31:40))

# Define the column(s) you do not want to impute
# Create a default method vector using `mice` defaults and set `""` for excluded columns
methods <- make.method(vic_reptile_traits_to_impute, defaultMethod = c("rf", "rf", "polyreg", "polr"))

pred.mat = as.matrix(c(0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
pred.mat = pred.mat[1:30,1]
pred.mat = matrix(rep(pred.mat, each = length(pred.mat)), nrow = length(pred.mat), byrow = TRUE)

imputed_data <- mice(vic_reptile_traits_to_impute, 
                     m = 5, # How many imputations
                     predictorMatrix = pred.mat,
                     method = methods, # Method used (needs to be specific to continuous vs. categorical)
                     maxit = 50, # Number of iterations for algorithm 
                     seed = 123) # Set Seed


# Bring imputed data back into the dataset 
vic_reptile_traits_imputed <- complete(imputed_data)

vic_reptile_traits_imputed <- vic_reptile_traits %>%
  select(-all_of(names(vic_reptile_traits_to_impute)[2:30])) %>%  # Remove columns in df1 that need replacement
  left_join(vic_reptile_traits_imputed, by = "Taxon_ID")  %>%         # Join with df2 to bring in the replacement columns
  select(all_of(names(vic_fauna_traits)))

# Save
saveRDS(vic_reptile_traits_imputed, "data_clean/vic_reptile_traits_imputed.Rds")


################# 
#### Frogs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################
vic_frog_traits = vic_fauna_traits %>% filter(Taxa_Group=="Amphibians") %>% droplevels()

#vic_frog_traits = readRDS("data_clean/vic_frog_traits.Rds")

#mice::md.pattern(vic_frog_traits[,1:17])
#pred.mat = as.matrix(c(0,0,0,0,1,1,0,0,0,1,1,1,1,1,1,1,1))
#pred.mat = matrix(rep(pred.mat, each = length(pred.mat)), nrow = length(pred.mat), byrow = TRUE)

# Pick the traits we actually have enough data to impute - the rest will stay NAs
vic_frog_traits_to_impute = vic_frog_traits %>% select(c(1:8, Mass_g, n_offspring_year))

pred.mat = as.matrix(c(0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
pred.mat = pred.mat[1:ncol(vic_frog_traits_to_impute),1]
pred.mat = matrix(rep(pred.mat, each = length(pred.mat)), nrow = length(pred.mat), byrow = TRUE)

# Define the column(s) you do not want to impute
# Create a default method vector using `mice` defaults and set `""` for excluded columns
methods <- make.method(vic_frog_traits_to_impute, defaultMethod = c("rf", "rf", "polyreg", "polr"))

imputed_data <- mice(vic_frog_traits_to_impute, 
                     m = 5, # How many imputations
                     predictorMatrix = pred.mat,
                     method = methods, # Method used (needs to be specific to continuous vs. categorical)
                     maxit = 50, # Number of iterations for algorithm 
                     seed = 123) # Set Seed


# Bring imputed data back into the dataset 
vic_frog_traits_imputed <- complete(imputed_data)
#vic_frog_traits_imputed <- cbind(vic_frog_traits_imputed, vic_frog_traits[,44:52])


vic_frog_traits_imputed <- vic_frog_traits %>%
  select(-all_of(names(vic_frog_traits_to_impute)[2:10])) %>%  # Remove columns in df1 that need replacement
  left_join(vic_frog_traits_imputed, by = "Taxon_ID")  %>%         # Join with df2 to bring in the replacement columns
  select(all_of(names(vic_fauna_traits)))


# Save
saveRDS(vic_frog_traits_imputed, "data_clean/vic_frog_traits_imputed.Rds")

# Join back together
vic_fauna_traits_imputed = rbind(vic_mammal_traits_imputed, vic_bird_traits_imputed, 
                                 vic_reptile_traits_imputed, vic_frog_traits)

##################################
### Save
saveRDS(vic_fauna_traits_imputed, "data_clean/vic_fauna_traits_imputed.Rds")
trait.names = data.frame(traits = names(vic_fauna_traits_imputed))
write.csv(trait.names, "data_clean/trait.names.csv")
