# ------------------------------------------------------------------------------
# Script Name: Step 1_Get species lists
# Author: Billy Geary
# Date Created: 17/12/2024
# Last Modified: 17/12/2024
# Purpose: This script creates an initial species list to act as the source of truth for the EMV Project
#          analysing species vulnerability to fire based on trait data.
#          The script reads in the Species data that DEECA has HDMs for and corrects any taxonomy issues based on the
#          FFG Threatened List as of June 2024. 
# Outputs: A final species list to be used in the EMV project, along with key taxonomic and threatened status information. 
# ------------------------------------------------------------------------------

## Step 1: Get fauna species
library(tidyverse)
library(readxl)

# Species list
species <- readxl::read_xlsx("data_clean/SMP_inputfiles_v3.0_final.xlsx", sheet="species")
thr.list <- readxl::read_xlsx("data_raw/DEECA data/FFG_Threatened_List_June_2024.xlsx")


# Pick the columns & filter to fauna
species <- species %>% select(TAXON_ID, SCIENTIFIC_NAME, COMMON_NAME, TaxaGroup, TAXON_TYPE, FFG_ACT_STATUS)
fauna <- species %>% filter(TaxaGroup %in% c("Mammals", "Amphibians", "Birds", "Reptiles"))
thr.list = thr.list %>% select(`Taxon ID`, `Scientific Name`, `Extinction Risk`, `Category of Threat`)
fauna <- left_join(fauna, thr.list, by=c("TAXON_ID"="Taxon ID"))

fauna = fauna %>%
  mutate(Sci_Name = ifelse(is.na(`Scientific Name`), SCIENTIFIC_NAME, `Scientific Name`),
         Extinction_Risk = ifelse(is.na(`Extinction Risk`), "NotListed", `Extinction Risk`),
         FFG_Status = ifelse(is.na(`Category of Threat`), "NotListed", `Category of Threat`)) %>%
  select(TAXON_ID, Sci_Name, COMMON_NAME, TaxaGroup, TAXON_TYPE, Extinction_Risk, FFG_Status)

# Fix up the taxonomy for species

mammal.fixes = data.frame(Sci_Name = c("Mormopterus sp. 4", "Mormopterus sp. 2", "Petaurus australis australis",
                                       "Mormopterus sp. 3", "Cercartetus concinnus minor", "Sminthopsis murina murina",
                                       "Tadarida australis", "Isoodon obesulus obesulus", "Miniopterus orianae oceanensis",
                                       "Potorous tridactylus trisulcatus", "Nyctophilus corbeni", "Dasyurus maculatus maculatus",
                                       "Rhinolophus megaphyllus megaphyllus", "Miniopterus orianae bassanii", "Osphranter robustus robustus",
                                       "Mastacomys fuscus mordicus", "Antechinus minimus maritimus", "Notomys mitchelli"),
                          Sci_Name_Alt = c("Mormopterus planiceps", "Mormopterus norfolkensis", "Petaurus australis",
                                           "Ozimops petersi", "Cercartetus concinnus", "Sminthopsis murina",
                                           "Austronomus australis", "Isoodon obesulus", "Miniopterus schreibersii",
                                           "Potorous tridactylus", "Nyctophilus corbeni", "Dasyurus maculatus",
                                           "Rhinolophus megaphyllus", "Miniopterus schreibersii", "Macropus robustus",
                                           "Mastacomys fuscus", "Antechinus minimus", "Notomys mitchellii"))

bird.fixes = data.frame(Sci_Name = c("Lichenostomus chrysops", "Lichenostomus penicillatus", "Lichenostomus leucotis", 
                                     "Lichenostomus virescens", "Lichenostomus ornatus", "Lichenostomus fuscus", 
                                     "Artamus leucorynchus", "Ardea ibis", "Phylidonyris albifrons", 
                                     "Phylidonyris melanops", "Phylidonyris pyrrhoptera", "Coracina tenuirostris", 
                                     "Cracticus tibicen", "Chrysococcyx basalis", "Cacomantis pallidus", 
                                     "Tyto javanica", "Chrysococcyx lucidus", "Accipiter cirrhocephalus", 
                                     "Climacteris picumnus victoriae", "Turnix varia", "Coturnix ypsilophora australis", 
                                     "Calyptorhynchus funereus", "Chrysococcyx osculans", "Porzana tabuensis", 
                                     "Pardalotus punctatus punctatus", "Threskiornis molucca", "Psephotus varius", 
                                     "Todiramphus pyrropygia pyrropygia", "Sericornis magnirostris", 
                                     "Cinclosoma castanotus", "Chlidonias hybridus javanicus", 
                                     "Larus pacificus pacificus", "Eolophus roseicapillus", "Cormobates leucophaeus", 
                                     "Cheramoeca leucosternus", "Megalurus gramineus", "Chroicocephalus novaehollandiae", 
                                     "Nycticorax caledonicus hillii", "Sugamel niger", "Gallirallus philippensis", 
                                     "Alcedo azurea", "Philomachus pugnax", "Amytornis striatus striatus", 
                                     "Pandion cristatus", "Butorides striatus", "Porzana pusilla palustris", 
                                     "Ardea alba modesta", "Ardea intermedia plumifera", "Lophocroa leadbeateri", 
                                     "Egretta garzetta nigripes", "Acanthiza iredalei hedleyi", 
                                     "Gelochelidon nilotica macrotarsa", "Polytelis anthopeplus monarchoides", 
                                     "Calyptorhynchus banksii graptogyne", "Lichenostomus melanops cassidix", 
                                     "Dupetor flavicollis"),
                        Sci_Name_Alt = c("Caligavis chrysops", "Ptilotula penicillata", "Nesoptilotis leucotis", 
                                         "Gavicalis virescens", "Ptilotula ornata", "Ptilotula fusca", 
                                         "Artamus leucoryn", "Bubulcus ibis", "Purnella albifrons", 
                                         "Gliciphila melanops", "Phylidonyris pyrrhopterus", "Edolisoma tenuirostre", 
                                         "Gymnorhina tibicen", "Chalcites basalis", "Heteroscenes pallidus", 
                                         "Tyto alba", "Chalcites lucidus", "Accipiter cirrocephalus", #####
                                         "Climacteris picumnus", "Turnix varius", "Synoicus ypsilophorus", 
                                         "Zanda funerea", "Chalcites osculans", "Zapornia tabuensis", 
                                         "Pardalotus punctatus", "Threskiornis moluccus", "Psephotellus varius", 
                                         "Todiramphus pyrrhopygius", "Sericornis magnirostra", #####
                                         "Cinclosoma castanotum", "Chlidonias hybrida", 
                                         "Larus pacificus", "Eolophus roseicapilla", "Cormobates leucophaea", 
                                         "Cheramoeca leucosterna", "Poodytes gramineus", "Larus novaehollandiae", 
                                         "Nycticorax caledonicus", "Sugomel nigrum", "Hypotaenidia philippensis", ####
                                         "Ceyx azureus", "Calidris pugnax", "Amytornis striatus", 
                                         "Pandion haliaetus", "Butorides striata", "Zapornia pusilla", 
                                         "Ardea alba", "Ardea intermedia", "Cacatua leadbeateri", 
                                         "Egretta garzetta", "Acanthiza iredalei", 
                                         "Gelochelidon nilotica", "Polytelis anthopeplus", 
                                         "Calyptorhynchus banksii", "Lichenostomus melanops", 
                                         "Ixobrychus flavicollis"))


reptile.fixes = data.frame(Sci_Name = c("Pseudemoia form cryodoma/pagenstecheri", "Eulamprus tympanum tympanum", 
                                        "Liopholis whitii GROUP", "Ramphotyphlops bituberculatus", "Parasuta nigriceps", 
                                        "Parasuta flagellum", "Ramphotyphlops bicolor", "Ramphotyphlops nigrescens", 
                                        "Parasuta dwyeri", "Ramphotyphlops proximus", 
                                        "Pseudemoia pagenstecheri (Volcanic Plains)", "Egernia saxatilis intermedia", 
                                        "Diplodactylus damaeus", "Amphibolurus nobbi", "Rankinia diemensis (Grampians)", 
                                        "Acritoscincus duperreyi", "Niveoscincus coventryi", "Niveoscincus metallicus", 
                                        "Acritoscincus platynotus", "Anepischtos maccoyi", "Rhinoplocephalus nigrescens", 
                                        "Intellagama lesueurii howittii", "Morelia spilota metcalfei", 
                                        "Eulamprus tympanum marnieae", "Parasuta spectabilis", "Morelia spilota spilota", 
                                        "Rankinia diemensis (Anglesea)"),
                           Sci_Name_Alt = c("Pseudemoia cryodoma", "Eulamprus tympanum", 
                                            "Liopholis whitii", "Anilios bituberculatus", "Suta nigriceps", 
                                            "Suta flagellum", "Anilios bicolor", "Anilios nigrescens", 
                                            "Suta dwyeri", "Anilios proximus", 
                                            "Pseudemoia pagenstecheri", "Egernia saxatilis", 
                                            "Lucasium damaeum", "Diporiphora nobbi", "Rankinia diemensis", 
                                            "Acritoscincus duperreyi", "Carinascincus coventryi", "Carinascincus metallicus", 
                                            "Acritoscincus platynotus", "Anepischetosia maccoyi", "Cryptophis nigrescens", 
                                            "Intellagama lesueurii", "Morelia spilota", 
                                            "Eulamprus tympanum", "Suta spectabilis", "Morelia spilota", 
                                            "Rankinia diemensis"))

frog.fixes = data.frame(Sci_Name = c("Litoria verreauxii verreauxii", "Litoria lesueuri", "Limnodynastes dumerilii dumerilii", 
                                     "Limnodynastes dumerilii insularis", "Limnodynastes dumerilii variegatus", 
                                     "Litoria watsoni", "Litoria verreauxii alpina", "Litoria booroolongensis"), 
                        Sci_Name_Alt =  c("Litoria verreauxii", "Litoria lesueurii", "Limnodynastes dumerilii", 
                                          "Limnodynastes dumerilii", "Limnodynastes dumerilii", 
                                          "Litoria littlejohni", "Litoria verreauxii", "Litoria booroolongensis"))

fixes = rbind(mammal.fixes, bird.fixes, reptile.fixes, frog.fixes)

fauna = left_join(fauna, fixes, by = "Sci_Name")
fauna$Sci_Name_Alt = ifelse(is.na(fauna$Sci_Name_Alt), fauna$Sci_Name, fauna$Sci_Name_Alt)
fauna$Sci_Name_Alt = ifelse(fauna$TAXON_ID==11136, "Petaurus australis", fauna$Sci_Name_Alt)


# Save 
write.csv(fauna, "data_clean/species_list_fauna.csv")
