
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

### Data Preparation for ToxPi score ###

library(readxl)
library(dplyr)

HBM <- read_excel('1-Overlap_EU_US_CA_N=76.xlsx',
                  sheet = 1)
### LD50 Prediction data ####
load('8.1-HBM_LD50_Pred_N=47.RData') #47 compounds
LD50_Pred <- LD50_Pred %>%
  select(-PREFERRED_NAME)
unique(LD50_Pred$DTXSID)

### Chemical Properties data ####
load('8.2-HBM_Kow_Koa_Pred_N=74.RData')

# Filter data, OPERA 1st option, EPISUITE 2nd option
# Selection Rules
Kow_Koa <- Kow_Koa %>%
  group_by(DTXSID, NAME) %>%
  # Use arrange to ensure "opera" comes first, then "episuite" comes second
  arrange(factor(SOURCE, levels = c('OPERA 2.6', 'EPISUITE')), 
          .by_group = TRUE) %>%
  # Select the first record of each group
  slice(1) %>%
  ungroup()

unique(Kow_Koa$DTXSID)

Kow_Koa <- Kow_Koa %>%
  select(DTXSID, NAME, VALUE)

library(tidyr)
#Convert to wide format式
Kow_Koa_1 <- Kow_Koa %>%
  # 确保 value 列是数值型
  mutate(VALUE = as.numeric(VALUE)) %>%
  # Make sure the value column is numeric
  pivot_wider(
    names_from = NAME,  # Create a new column from the name column
    values_from = VALUE  # Populate these new columns from the value column
  ) #74 compounds

### CTV Pre+Exp data ####
load('8.3-HBM_CTV_Merge.RData')

# Converting a data frame
CTV_1 <- CTV %>%
  select(-Unit, -Type) %>%
  pivot_wider(names_from = Tox_parameter, values_from = Value) %>%
  select(-PREFERRED_NAME) #63 compounds

CTV_1 <- CTV_1 %>%
  mutate(across(3:10, ~ as.numeric(.)))

### Cyto_AC50 ####
load('8.5-HBM_Cyto_AC50_Min_N=45.RData')

AC50_Min <- HBM_Cyto_AC50_Min %>%
  select(-ac50_2731, -ac50_3070, -ac50_3098, -ac50_1971, 
         -ac50_3163, -ac50_2031) %>%
  rename(DTXSID = dtxsid) #45 compounds

### Cyto_%Activity Assay ####
load('8.6-HBM_Cyto_%ActiveAssay_N=45.RData')

HIT_Assay_Ratio <- HBM_Cyto_HIT_Assay_Ratio %>%
  select(dtxsid, ratio) %>%
  rename(HIT_Ratio = ratio,
         DTXSID = dtxsid) #45 compounds


#### Merge all data ####
HBM_1 <- HBM %>%
  select(PREFERRED_NAME, DTXSID, CASRN, InChIKey, Chemical_Group)

HBM_Tox <- HBM_1 %>%
  left_join(Kow_Koa_1, by = 'DTXSID') %>%
  left_join(LD50_Pred, by = 'DTXSID') %>%
  left_join(CTV_1, by = 'DTXSID') %>%
  left_join(HIT_Assay_Ratio, by = 'DTXSID') %>%
  left_join(AC50_Min, by = 'DTXSID') 

save(HBM_Tox, file = '9-HBM_Tox_N=76.RData')









