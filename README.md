## Read Me

Repository for scripts to synthesise trait data related to fire sensitivity for mammal, bird, reptile and amphibian species in Victoria, Australia. 

### Versions

- Version 1.0: Current commit compiles traits for all species and constructs species clusters based on 1) SMP Benefit of Reducing Fire Freuncy and 2) FAME Fire Response curves. There remains some significant gaps in the frog and reptile datasets that would be good to continue to fill. Clustering remains a work in progress, but trait datasets are reasonably complete based on what's available publicly (though still gappy).

- Version 0.1: Initial commit for trait database compilation and initial species clustering based on SMP EE data. More to be done in terms of trait compilation and testing of alternative fire response measures (e.g. FAME). Also need to fix gaps in reptile and frog traits (e.g. diet, nesting, stratum) before results are sensible with trait importance scores, and also find better method for small sample size (e.g. frogs). 

### Scripts

- Step 1_Get species list

- Step 2_Biogeographic calculations

- Step 3_Compile FAME data

- Step 4_Compile trait data

- Step 5_Trait imputation

- Step 6_Identification of important traits

- Step 7_Map vulnerability spatially