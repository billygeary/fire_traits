library(tidyverse)
library(lme4)

# SMP EE Data 
# Fire Frequency
freq.fire = readRDS("data_raw/DEECA data/SMP_EE_planBurn.Rds") # DEECA EE Round 2
ee.data = freq.fire %>% 
  filter(action == "FrequentFRB" & tempscale == 50) %>%
  select(TaxonCode, group, expert, round, meanBen) %>%
  filter(group!="Plants") %>% mutate(TaxonCode = as.factor(TaxonCode))

ee.model = lmer(meanBen ~ TaxonCode + group + (1|expert) + (1|round), data = ee.data)
new.dat = ee.data %>% select(TaxonCode, group) %>% distinct() %>% drop_na()

ee.predict = data.frame(meanben_firefreq = predict(ee.model, newdata = new.dat, re.form = NA))

smp_ben_data = cbind(new.dat, ee.predict)


# Severe Fires
sev.fire = readRDS("data_raw/DEECA data/SMP_EE_futureFire.Rds")
# sev.fire = sev.fire %>% 
#   select(species, TaxonCode, scenarioID, lat, long, refmean, FutureFire, expert, group) %>% 
#   distinct() %>%
#   filter(group!="Plants") %>% 
#   pivot_wider(id_cols=c(species, TaxonCode, lat, long, expert, group), 
#               names_from = FutureFire, values_from = refmean, values_fn = mean) %>%
#   mutate(FutureImpact = `0` - `3`)# %>% drop_na(FutureImpact) 
# 
# 

sev.fire = sev.fire %>% 
  group_by(species, TaxonCode, FutureFire) %>% 
  summarise(meanDN = mean(refmean)) %>%
  pivot_wider(id_cols=c(species, TaxonCode), 
              names_from = FutureFire, values_from = meanDN) %>%
  mutate(FutureImpact = `0` - `3`) %>% drop_na(FutureImpact) %>% ungroup() %>%
  select(TaxonCode, FutureImpact)  %>% mutate(TaxonCode = as.factor(TaxonCode))

### Join Together

smp_ben_data = left_join(smp_ben_data, sev.fire, by = "TaxonCode")

write.csv(smp_ben_data, "data_clean/mean_smp_fire_data.csv")
