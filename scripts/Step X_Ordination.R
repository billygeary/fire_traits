## Trait Syndrome Analysis

## Start with an ordination of mammal traits
vic_mammal_traits_imputed = readRDS("data_clean/vic_mammal_traits_imputed.Rds")
dat = vic_mammal_traits_imputed

dat$terrestrial_volant = as.numeric(dat$terrestrial_volant)
dat = select(dat, c("Mass_g", 
                    "home_range_km2", 
                    "dispersal_km", 
                    "habitat_breadth_n",
                    #"terrestrial_volant", 
                    "hibernation_torpor",
                    "fossoriality", 
                    "det_diet_breadth_n", 
                    #"dphy_invertebrate", "dphy_vertebrate", "dphy_plant",  
                    "max_longevity_d", 
                    "n_offspring_year", 
                    #"Aerial", 
                    "Arboreal_Insessorial", "Terrestrial",
                    "biggest_patch_size", "n_patches", "nned",
                    "pyrome_breadth", "dominant_pyrome",
                    "meanBen_FireFreq")) %>%
  mutate(Mass_g = Mass_g/1000) %>% drop_na()

### Correlation Plot to visualise ###
cors = cor(dat)
corrplot::corrplot(cors)

## Get ready for ordination
dat = select(dat, c("Mass_g", 
                    #"home_range_km2", 
                    #"dispersal_km", 
                    #"habitat_breadth_n",
                    #"terrestrial_volant", 
                    #"hibernation_torpor",
                    #"fossoriality", 
                    #"det_diet_breadth_n", 
                    #"dphy_invertebrate", "dphy_vertebrate", "dphy_plant",  
                    #"max_longevity_d", 
                    "n_offspring_year", "pyrome_breadth",
                    "Aerial", "Arboreal_Insessorial", "Terrestrial",
                    "meanBen_FireFreq",  "FutureImpact")) %>%
  #mutate(Mass_g = Mass_g/1000) %>% 
  drop_na()

gd <- FD::gowdis(dat)

# perform principal coordinates analysis (PCoA)
diet_pco <- ade4::dudi.pco(gd, scannf = FALSE)

pc_diet <- diet_pco$tab
summary(diet_pco)

# principle component axes
pcomp_diet <- as.data.frame(diet_pco$tab[,1:2]) %>% 
  tibble::rownames_to_column(var = "species")
names(pcomp_diet) <- c("species", "A1", "A2")

# diet category projection
n <- nrow(dat)
points_stand <- scale(diet_pco$tab[,1:2])
S <- cov(dat, points_stand)
U <- S %*% diag((diet_pco$eig[1:2]/(n-1))^(-0.5))
colnames(U) <- colnames(diet_pco$tab[,1:2])

U <- as.data.frame(U) %>% 
  mutate(trait = names(dat))

# scale diet category arrows
mult <- min(
  (max(pcomp_diet$A2) - min(pcomp_diet$A2)/(max(U$A2)-min(U$A2))),
  (max(pcomp_diet$A1) - min(pcomp_diet$A1)/(max(U$A1)-min(U$A1)))
)

U <- U %>% 
  mutate(v1 = 0.03 * mult * A1) %>% 
  mutate(v2 = 0.03 * mult * A2)

# plot diet PCoA
pcoa_diet <- ggplot(pcomp_diet, aes(A1, A2)) +
  # set up plot
  geom_hline(yintercept = 0, size = 0.2, lty = 2, colour = "grey") + 
  geom_vline(xintercept = 0, size = 0.2, lty = 2, colour = "grey") +
  # add origin lines
  #geom_text(alpha = 0.4, size = 1, aes(label = binomial))
  geom_point() +
  # add species
  coord_equal() +
  geom_segment(data = U, aes(x = 0, y = 0, xend = v1, yend = v2), arrow = arrow(length = unit(0.2, "cm")), colour = "darkgrey") +
  # add arrows
  geom_text(data = U, aes(x = v1, y = v2, label = trait), size = 4, colour = "darkgrey",
            nudge_y = c(rep(0, 6), 0.005, 0.005, 0.0005, -0.004), 
            nudge_x = c(0.005, rep(0, 7), -0.009, 0)) +
  # add arrow labels
  labs(x = "PC1 (52.62%)", y = "PC2 (28.41%)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 10))

pcoa_diet


m <- lm(meanBen_FireFreq~ scale(pyrome_breadth) + scale(Mass_g) + 
                           scale(n_offspring_year), data=dat)
summary(m)

## Start with an ordination of bird traits
vic_bird_traits_imputed = readRDS("data_clean/vic_bird_traits_imputed.Rds")
dat = vic_bird_traits_imputed

### Correlation Plot to visualise ###
dat = dat %>% drop_na(meanBen_FireFreq)
cors = cor(dat[,14:38])
corrplot::corrplot(cors)

dat = select(dat, c("Mass_g", 
                    #"home_range_km2", 
                    #"dispersal_km", 
                    #"habitat_breadth_n",
                    #"terrestrial_volant", 
                    #"hibernation_torpor",
                    #"fossoriality", 
                    #"det_diet_breadth_n", 
                    #"dphy_invertebrate", "dphy_vertebrate", "dphy_plant",  
                    #"max_longevity_d", 
                    "litters_or_clutches_per_y", 
                    "pyrome_breadth",
                    "Aerial", "Arboreal_Insessorial", "Terrestrial",
                    "Generalist", "Aquatic",
                    "meanBen_FireFreq","FutureImpact")) %>%
  mutate(Mass_g = Mass_g/1000) %>% drop_na()

gd <- FD::gowdis(dat)

# perform principal coordinates analysis (PCoA)
diet_pco <- ade4::dudi.pco(gd, scannf = FALSE)

pc_diet <- diet_pco$tab
summary(diet_pco)

# principle component axes
pcomp_diet <- as.data.frame(diet_pco$tab[,1:2]) %>% 
  tibble::rownames_to_column(var = "species")
names(pcomp_diet) <- c("species", "A1", "A2")

# diet category projection
n <- nrow(dat)
points_stand <- scale(diet_pco$tab[,1:2])
S <- cov(dat, points_stand)
U <- S %*% diag((diet_pco$eig[1:2]/(n-1))^(-0.5))
colnames(U) <- colnames(diet_pco$tab[,1:2])

U <- as.data.frame(U) %>% 
  mutate(trait = names(dat))

# scale diet category arrows
mult <- min(
  (max(pcomp_diet$A2) - min(pcomp_diet$A2)/(max(U$A2)-min(U$A2))),
  (max(pcomp_diet$A1) - min(pcomp_diet$A1)/(max(U$A1)-min(U$A1)))
)

U <- U %>% 
  mutate(v1 = 0.03 * mult * A1) %>% 
  mutate(v2 = 0.03 * mult * A2)

# plot diet PCoA
pcoa_diet <- ggplot(pcomp_diet, aes(A1, A2)) +
  # set up plot
  geom_hline(yintercept = 0, size = 0.2, lty = 2, colour = "grey") + 
  geom_vline(xintercept = 0, size = 0.2, lty = 2, colour = "grey") +
  # add origin lines
  #geom_text(alpha = 0.4, size = 1, aes(label = binomial))
  geom_point() +
  # add species
  coord_equal() +
  geom_segment(data = U, aes(x = 0, y = 0, xend = v1, yend = v2), arrow = arrow(length = unit(0.2, "cm")), colour = "darkgrey") +
  # add arrows
  geom_text(data = U, aes(x = v1, y = v2, label = trait), size = 4, colour = "darkgrey",
            nudge_y = c(rep(0, 6), 0.005, 0.005, 0.0005, -0.004), 
            nudge_x = c(0.005, rep(0, 7), -0.009, 0)) +
  # add arrow labels
  labs(x = "PC1 (42.84%)", y = "PC2 (18.14%)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 10))

pcoa_diet
