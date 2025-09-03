
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

##### 1. ChemPro and LD50 data for HBM 76 compounds ####
### Download Physicochemical Property Values, TEST Model Predictions, OPERA Model Predictions data from CompTox
# Download data: Liver_BatchSearch_CompTox.xlsx

library(readxl)
### Extract 'ORAL_RAT_LD50' from TEST ####
TEST <- read_excel('Liver_BatchSearch_CompTox.xlsx',
                   sheet = 'Main Data')
unique(TEST$DTXSID) #94 compounds

LD50_Pred <- TEST %>%
  select(DTXSID, PREFERRED_NAME, `ORAL_RAT_LD50_MOL/KG_TEST_PRED`) %>%
  filter(!is.na(`ORAL_RAT_LD50_MOL/KG_TEST_PRED`)) #79 compounds

save(LD50_Pred, file = '6.1-Liver_LD50_Pred_N=79.RData')

### Extract 'LogKow', 'LogKoa' from ChemPro ####
ChemPro <- read_excel('Liver_BatchSearch_CompTox.xlsx',
                      sheet = 'Chemical Properties')
unique(ChemPro$DTXSID) #90 compounds

Kow_Koa <- ChemPro %>%
  select(-DTXCID) %>%
  filter(NAME %in% c('LogKoa: Octanol-Air', 'LogKow: Octanol-Water')) %>%
  filter(SOURCE %in% c('EPISUITE','OPERA 2.6'))

unique(Kow_Koa$DTXSID)
save(Kow_Koa, file = '6.2-Liver_Kow_Koa_Pred_N=88.RData')

##### 2. Conditional Toxicity Value (CTV) #####
### Download data from CTV website: https://toxvalue.org/6-CTV/Cover.php
### Raw data clean process in the HBM_CTV_RawData_Clean file
load('Liver_CTV_RawData_Clean/Liver_CTV_Pre_Exp.RData')

# Revise variable name
library(dplyr)
unique(CTV_Exp$Endpoint)
CTV_Exp_1 <- CTV_Exp %>%
  mutate(Tox_parameter = case_when(Endpoint == 'Reference Dose' ~ 'RfD',
                                   Endpoint == 'Reference Concentration'  ~ 'RfC',
                                   Endpoint == "Oral Slope Factor"  ~ 'OSF',
                                   Endpoint == "Cancer Potency Value"  ~ 'CPV',
                                   Endpoint == "Inhalation Unit Risk"  ~ 'IUR')
  ) %>%
  mutate(Type = 'Experiment') %>%
  rename( Value = `Toxicity value`) %>%
  select(-Endpoint)

unique(CTV_Prediction$Model.Name)
CTV_Pred_1 <- CTV_Prediction %>%
  mutate(Tox_parameter = case_when(Model.Name == "CTV Reference Dose (RfD)" ~ 'RfD',
                                   Model.Name == "CTV Reference Dose NO(A)EL" ~ 'NOAEL',
                                   Model.Name == "CTV Reference Dose (RfD) BMD" ~ 'BMD',
                                   Model.Name == "CTV Reference Dose (RfD) BMDL" ~ 'BMDL',
                                   Model.Name == "CTV Reference Concentration (RfC)" ~ 'RfC',
                                   Model.Name == "CTV Oral Slope Factor (OSF)" ~ 'OSF',
                                   Model.Name == "CTV Cancer Potency Value (CPV)" ~ 'CPV',
                                   Model.Name == "CTV Inhalation Unit Risk (IUR)" ~ 'IUR')
  ) %>%
  mutate(Type = 'Prediction') %>%
  dplyr::rename(Value = Prediction) %>%
  select(-Model.Name)

# Merge database Exp+Pred
CTV <- rbind(CTV_Exp_1, CTV_Pred_1)

save(CTV, file = '6.3-Liver_CTV_Merge.RData')

##### 3. CompTox Cytotoxicity Data #####

### Select Cytotoxicity Data
load('3-Human_Liver_Cytotoxicity_bio_merge.RData')

## Load HBM 76 compounds
library(readxl)
Liver_CS <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                     sheet = 1)
# Extract DTXID and save as Values
Liver_DTXID <- Liver_CS %>%
  select(DTXSID) %>%
  distinct() %>%
  pull(DTXSID) #94 DTXSID

## Select HBM_compounds using dtxsid
Liver_cyto_bio <- cyto_bio_merge %>%
  filter(dtxsid %in% Liver_DTXID) # 2649 records

dtxsid_list <- Liver_cyto_bio$dtxsid[1:2649] # 82 compounds
unique(dtxsid_list) # unique DTXSID 

# Save data
save(Liver_cyto_bio, file = '6.4-Liver_Cytotoxicity.RData')

### Clean Cytotoxicity Data
# Load Liver_cytotoxicity data
load('6.4-Liver_Cytotoxicity.RData')
Liver_Cyto <- Liver_cyto_bio #2649 records

# Extract variables from data
Liver_Cyto_1 <- subset(Liver_Cyto, 
                     select = c(ac50, hitcall, dtxsid, maxMeanConc, aeid, spid))

# Add a new variable
library(dplyr)
Liver_Cyto_2 <- Liver_Cyto_1 %>%
  mutate(HIT_CALL = factor(case_when(
    hitcall < 0.9 ~ "Inactive",
    hitcall >= 0.9 ~ "Active"
  )))

## Output Cytotoxicity AC50 value

# Extract AEID, AC50, DTXSID
Liver_Cyto_AC50 <- subset(Liver_Cyto_2,
                        select = c(dtxsid, aeid, spid, ac50)) #2649 records

# Use the pivot_wider function to change the aeid value into a column name
library(tidyr)
library(dplyr)

# 
Liver_Cyto_AC50_unique <- Liver_Cyto_AC50 %>%
  distinct(dtxsid, aeid, ac50, .keep_all = TRUE) #2407 records

# 
Liver_Cyto_AC50_unique <- Liver_Cyto_AC50_unique %>%
  group_by(dtxsid, aeid, spid) %>%
  mutate(record_id = row_number()) %>%
  ungroup()

# 
Liver_Cyto_AC50_Assay <- Liver_Cyto_AC50_unique %>%
  pivot_wider(
    names_from = aeid,
    values_from = ac50,
    names_prefix = "ac50_"
  ) %>%
  arrange(dtxsid, record_id) %>%
  select(-record_id) #376 records, 99 cytotoxicity assays

unique(Liver_Cyto_AC50_Assay$dtxsid) #82 DTXSID

### Select the min_AC50 for each compound across different 'spid'#### 
df <- Liver_Cyto_AC50_Assay %>%
  select(-spid) 

# Convert empty strings ("") to NA
library(dplyr)
df <- df %>%
  mutate(across(where(is.character), ~ na_if(., "")))

# Group by dtxsid and select the smallest value in each group, ignoring NA
df <- df %>%
  group_by(dtxsid) %>%
  summarize(across(everything(), ~ {
    if (all(is.na(.))) {
      NA  
    } else {
      min_value <- min(., na.rm = TRUE)
      if (is.infinite(min_value)) NA else min_value
    }
  }))

Liver_Cyto_AC50_Min <- df
save(Liver_Cyto_AC50_Min, file = '6.5-Liver_Cyto_AC50_Min_N=82.RData')

#### Output Cytotoxicity HIT_CALL, calculate %Active assays ####

# Extract aeid, dtxsid, HIT_CALL
Liver_Cyto_HIT <- subset(Liver_Cyto_2,
                       select = c(dtxsid, aeid, spid, HIT_CALL)) #2649 records

# Remove duplicate values
Liver_Cyto_HIT_unique <- Liver_Cyto_HIT %>%
  distinct(dtxsid, aeid, HIT_CALL, .keep_all = TRUE) #1679 records

# HIT_CALL has duplicate values, so you need to create a unique identifier for each dtxsid and aeid, spid combination
Liver_Cyto_HIT_unique <- Liver_Cyto_HIT_unique %>%
  group_by(dtxsid, aeid, spid) %>%
  mutate(record_id = row_number()) %>%
  ungroup()

# Use pivot_wider and expand the HIT_CALL values into separate records
Liver_Cyto_HIT_Assay <- Liver_Cyto_HIT_unique %>%
  pivot_wider(
    names_from = aeid,
    values_from = HIT_CALL,
    names_prefix = "HIT_"
  ) %>%
  arrange(dtxsid, record_id) %>%
  select(-record_id) #179 records

unique(Liver_Cyto_HIT_Assay$dtxsid) #82 DTXSID

### Calculate the percentage of active assays

# Convert factor variables to numeric form (active is 1, inactive is 0)
Liver_Cyto_HIT_Assay_numeric <- Liver_Cyto_HIT_Assay %>%
  mutate(across(starts_with("HIT"), ~ ifelse(. == "Active", 1, 0)))

# 
str(Liver_Cyto_HIT_Assay_numeric)
# 
print(Liver_Cyto_HIT_Assay_numeric)


# Calculate the active percentage for each dtxsid
Liver_Cyto_HIT_Assay_percentage <- Liver_Cyto_HIT_Assay_numeric %>%
  group_by(dtxsid) %>%
  summarise(
    total_hits = sum(!is.na(across(starts_with("HIT_")))), # Total number of factor variables, ignoring NA
    active_count = sum(across(starts_with("HIT_"), ~ sum(. == 1, na.rm = TRUE))), # Calculate the sum of active values in all factor variables
    ratio = (active_count / total_hits) * 100, 
    .groups = 'drop' 
  )

# 
print(Liver_Cyto_HIT_Assay_percentage) # 82 records

# Merge calculation results, Delete repeated dtxsid record
Liver_Cyto_HIT_Assay_Ratio <- Liver_Cyto_HIT_Assay %>%
  left_join(Liver_Cyto_HIT_Assay_percentage, by = "dtxsid") %>%
  distinct(dtxsid, .keep_all = TRUE) #82 records

# 
print(Liver_Cyto_HIT_Assay_Ratio)

save(Liver_Cyto_HIT_Assay_Ratio, file = '6.6-Liver_Cyto_%ActiveAssay_N=82.RData')

#### Output Cytoxicity maxMeanConc #####
# Extract AEID, maxMeanConc, DTXSID
Liver_Cyto_Max <- subset(Liver_Cyto_2,
                       select = c(dtxsid, aeid, spid, maxMeanConc)) #2649 records

# 
library(tidyr)
library(dplyr)

# 
Liver_Cyto_Max_unique <- Liver_Cyto_Max %>%
  distinct(dtxsid, aeid, maxMeanConc, .keep_all = TRUE) #2490 records

# maxMeanConc has repeated values, so a unique identifier needs to be created for each dtxsid and aeid, spid combination
Liver_Cyto_Max_unique <- Liver_Cyto_Max_unique %>%
  group_by(dtxsid, aeid, spid) %>%
  mutate(record_id = row_number()) %>%
  ungroup()

# 
Liver_Cyto_Max_Assay <- Liver_Cyto_Max_unique %>%
  pivot_wider(
    names_from = aeid,
    values_from = maxMeanConc,
    names_prefix = "max_"
  ) %>%
  arrange(dtxsid, record_id) %>%
  select(-record_id) #219 records

unique(Liver_Cyto_Max_Assay$dtxsid) #82 DTXSID

library(writexl)
save(Liver_Cyto_Max_Assay, file = '6.7-Liver_Cyto_maxMeanConc_Assay.RData')



