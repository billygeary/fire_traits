# ------------------------------------------------------------------------------
# Script Name: Step 2_FAME Calculations
# Author: Billy Geary
# Date Created: 17/12/2024
# Last Modified: 17/12/2024
# Purpose: This script uses the DEECA Decision-support tool FAME's output data to estimate the response of each species to time since fire. 
#          Two datasets are used: 1) an empirical dataset derived from ERP 28_2 and 2) The FAME Future Fauna Occupancy expert opinion dataset. 
#          Slopes for each species are calculated. 
# Outputs: A set of csv files with the output of the calculations for each species. Species are ID'd using the VBA Taxon Codes.  
# ------------------------------------------------------------------------------

# Compile FAME data
library(MuMIn)
library(tidyverse)
fame = read.csv("data_raw/DEECA data/ERP28_2andNorthernMattleeBirdcurves/ERP28_2andNorthernMattleeBirdcurves.csv")
ids = na.omit(unique(fame$TAXON_ID))

data.out = data.frame()
for(s in ids){
  spp.data = filter(fame, TAXON_ID == s)
  mod = lm(Abund ~ YSF, data = spp.data)
  m = summary(mod)
  mod2 = lm(Abund ~ YSF+poly(YSF,2), data = spp.data)
  m2 = summary(mod2)
  
  # Find the peak / trough
  # Generate a sequence of YSF values from 0 to 300 (your range)
  ysf_seq <- seq(0, 400, length.out = 401)
  # Predict abundance using the model for each YSF value
  preds <- predict(mod2, newdata = data.frame(YSF = ysf_seq))
  # Find the YSF value that gives the maximum predicted abundance
  ysf_max <- ysf_seq[which.max(preds)]
  ysf_min <- ysf_seq[which.min(preds)]
  out = data.frame(TAXON_ID = s, 
                   lm_slope = m$coefficients[2,1],
                   lm_p_value = m$coefficients[2,4],
                   lm_aic = AICc(mod),
                   plm_slope = m2$coefficients[2,1],
                   plm_slope_p_value = m2$coefficients[2,4],
                   plm_poly = m2$coefficients[3,1],
                   plm_poly_pvalue = m2$coefficients[3,4],
                   plm_aic = AICc(mod2),
                   ysf_max = ysf_max,
                   ysf_min = ysf_min
                   )
  data.out = rbind(data.out, out)
}

data.out$model = ifelse(data.out$plm_aic < data.out$lm_aic, "poly", "lm")

data.out.cats = data.out %>% mutate(Sig = case_when(model == "poly" & 
                                              plm_slope_p_value < 0.05 | 
                                              plm_poly_pvalue < 0.05 ~ 1,
                                              model == "lm" & lm_p_value < 0.05 ~ 1,
                                              TRUE ~ 0)) %>%
  mutate(curve_shape = case_when(model == "lm" & Sig == 1 & lm_slope > 0 ~ "Positive",
                                 model == "lm" & Sig == 1 & lm_slope < 0 ~ "Negative",
                                 model == "poly" & Sig == 1 & plm_slope > 0 & plm_poly > 0 ~ "Positive_Trough",
                                 model == "poly" & Sig == 1 & plm_slope > 0 & plm_poly < 0 ~ "Positive_Peak",
                                 model == "poly" & Sig == 1 & plm_slope < 0 & plm_poly > 0 ~ "Negative_Trough",
                                 model == "poly" & Sig == 1 & plm_slope < 0 & plm_poly < 0 ~ "Negative_Peak",
                                 TRUE ~ "Neutral"))


fame.data = data.out.cats %>% select(TAXON_ID, curve_shape, lm_slope)
write.csv(fame.data, "data_clean/fame_erp_slopes.csv")

# Fame Expert Opinion Data

fame.ee = read.csv("data_raw/DEECA data/ERP28_2andNorthernMattleeBirdcurves/FFOResponsesLong.csv")

ids = na.omit(unique(fame.ee$TAXON_ID))
ee.data.out = data.frame()
for(s in ids){
  spp.data = filter(fame.ee, TAXON_ID == s)
  mod = lm(Abund ~ YSF, data = spp.data)
  m = summary(mod)
  mod2 = lm(Abund ~ YSF+poly(YSF,2), data = spp.data)
  m2 = summary(mod2)
  
  # Find the peak / trough
  # Generate a sequence of YSF values from 0 to 300 (your range)
  ysf_seq <- seq(0, 400, length.out = 401)
  # Predict abundance using the model for each YSF value
  preds <- predict(mod2, newdata = data.frame(YSF = ysf_seq))
  # Find the YSF value that gives the maximum predicted abundance
  ysf_max <- ysf_seq[which.max(preds)]
  ysf_min <- ysf_seq[which.min(preds)]
  out = data.frame(TAXON_ID = s, 
                   ee_lm_slope = m$coefficients[2,1],
                   ee_lm_p_value = m$coefficients[2,4],
                   ee_lm_aic = AICc(mod),
                   ee_plm_slope = m2$coefficients[2,1],
                   ee_plm_slope_p_value = m2$coefficients[2,4],
                   ee_plm_poly = m2$coefficients[3,1],
                   ee_plm_poly_pvalue = m2$coefficients[3,4],
                   ee_plm_aic = AICc(mod2),
                   ee_ysf_max = ysf_max,
                   ee_ysf_min = ysf_min
  )
  ee.data.out = rbind(ee.data.out, out)
}

ee.data.out$ee_model = ifelse(ee.data.out$ee_plm_aic < ee.data.out$ee_lm_aic, "poly", "lm")

ee.data.out.cats = ee.data.out %>% mutate(ee_Sig = case_when(ee_model == "poly" & 
                                                      ee_plm_slope_p_value < 0.05 | 
                                                      ee_plm_poly_pvalue < 0.05 ~ 1,
                                                    ee_model == "lm" & ee_lm_p_value < 0.05 ~ 1,
                                                    TRUE ~ 0)) %>%
  mutate(ee_curve_shape = case_when(ee_model == "lm" & ee_Sig == 1 & ee_lm_slope > 0 ~ "Positive",
                                 ee_model == "lm" & ee_Sig == 1 & ee_lm_slope < 0 ~ "Negative",
                                 ee_model == "poly" & ee_Sig == 1 & ee_plm_slope > 0 & ee_plm_poly > 0 ~ "Positive_Trough",
                                 ee_model == "poly" & ee_Sig == 1 & ee_plm_slope > 0 & ee_plm_poly < 0 ~ "Positive_Peak",
                                 ee_model == "poly" & ee_Sig == 1 & ee_plm_slope < 0 & ee_plm_poly > 0 ~ "Negative_Trough",
                                 ee_model == "poly" & ee_Sig == 1 & ee_plm_slope < 0 & ee_plm_poly < 0 ~ "Negative_Peak",
                                 TRUE ~ "Neutral"))


fame.data = ee.data.out.cats %>% select(TAXON_ID, ee_model, ee_lm_slope)

write.csv(fame.data, "data_clean/fame_ffo_slopes.csv")

### Combined data
### For each species, if we have erp data, use that, otherwise use ffos. 
erp = read.csv("data_clean/fame_erp_slopes.csv") %>% select(TAXON_ID, lm_slope)
ffo = read.csv("data_clean/fame_ffo_slopes.csv") %>% select(TAXON_ID, ee_lm_slope)

fame_slopes = full_join(erp, ffo)
fame_slopes$fame_lm_slope = ifelse(is.na(fame_slopes$lm_slope), fame_slopes$ee_lm_slope, fame_slopes$lm_slope)

write.csv(fame_slopes, "data_clean/fame_combined_slopes.csv")


