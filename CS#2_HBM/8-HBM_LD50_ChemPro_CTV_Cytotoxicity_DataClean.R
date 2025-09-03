
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

##### 1. ChemPro and LD50 data for HBM 76 compounds ####
### Download Physicochemical Property Values, TEST Model Predictions, OPERA Model Predictions data from CompTox
# Download data: HBM_BatchSearch_CompTox.xlsx

library(readxl)
### Extract 'ORAL_RAT_LD50' from TEST ####
TEST <- read_excel('HBM_BatchSearch_CompTox.xlsx',
                   sheet = 'Main Data')
unique(TEST$DTXSID) #76 compounds

LD50_Pred <- TEST %>%
  select(DTXSID, PREFERRED_NAME, `ORAL_RAT_LD50_MOL/KG_TEST_PRED`) %>%
  filter(!is.na(`ORAL_RAT_LD50_MOL/KG_TEST_PRED`)) #47 compounds

save(LD50_Pred, file = '8.1-HBM_LD50_Pred_N=47.RData')

### Extract 'LogKow', 'LogKoa' from ChemPro ####
ChemPro <- read_excel('HBM_BatchSearch_CompTox.xlsx',
                      sheet = 'Chemical Properties')
unique(ChemPro$DTXSID) #74 compounds

Kow_Koa <- ChemPro %>%
  select(-DTXCID) %>%
  filter(NAME %in% c('LogKoa: Octanol-Air', 'LogKow: Octanol-Water')) %>%
  filter(SOURCE %in% c('EPISUITE','OPERA 2.6'))

save(Kow_Koa, file = '8.2-HBM_Kow_Koa_Pred_N=74.RData')

##### 2. Conditional Toxicity Value (CTV) #####
### Download data from CTV website: https://toxvalue.org/6-CTV/Cover.php
### Raw data clean process in the HBM_CTV_RawData_Clean file
load('HBM_CTV_RawData_Clean/HBM_CTV_Pre_Exp.RData')

# Revise variable name
library(dplyr)
unique(CTV_Experiment$Endpoint)
CTV_Exp_1 <- CTV_Experiment %>%
  mutate(Tox_parameter = case_when(Endpoint == 'Reference Dose' ~ 'RfD',
                                   Endpoint == "Oral Slope Factor"  ~ 'OSF',
                                   Endpoint == "Cancer Potency Value"  ~ 'CPV',
                                   Endpoint == "Inhalation Unit Risk"  ~ 'IUR')
  ) %>%
  mutate(Type = 'Experiment') %>%
  dplyr::rename( Value = Toxicity_Value) %>%
  select(-Endpoint)

unique(CTV_Prediction$Model_Name)
CTV_Pred_1 <- CTV_Prediction %>%
  mutate(Tox_parameter = case_when(Model_Name == "CTV Reference Dose (RfD)" ~ 'RfD',
                                   Model_Name == "CTV Reference Dose NO(A)EL" ~ 'NOAEL',
                                   Model_Name == "CTV Reference Dose (RfD) BMD" ~ 'BMD',
                                   Model_Name == "CTV Reference Dose (RfD) BMDL" ~ 'BMDL',
                                   Model_Name == "CTV Reference Concentration (RfC)" ~ 'RfC',
                                   Model_Name == "CTV Oral Slope Factor (OSF)" ~ 'OSF',
                                   Model_Name == "CTV Cancer Potency Value (CPV)" ~ 'CPV',
                                   Model_Name == "CTV Inhalation Unit Risk (IUR)" ~ 'IUR')
  ) %>%
  mutate(Type = 'Prediction') %>%
  dplyr::rename(Value = Prediction) %>%
  select(-Model_Name)

# Merge database Exp+Pred
CTV <- rbind(CTV_Exp_1, CTV_Pred_1)

save(CTV, file = '8.3-HBM_CTV_Merge.RData')

##### 3. CompTox Cytotoxicity Data #####

### Select Cytotoxicity Data
load('4-Human_Cytotoxicity_bio_merge.RData')

## Load HBM 76 compounds
library(readxl)
HBM_CS <- read_excel('1-Overlap_EU_US_CA_N=76.xlsx',
                     sheet = 1)
# Extract DTXID and save as Values
HBM_DTXID <- HBM_CS %>%
  select(DTXSID) %>%
  distinct() %>%
  pull(DTXSID) #76 DTXSID

## Select HBM_compounds using dtxsid
HBM_cyto_bio <- cyto_bio_merge %>%
  filter(dtxsid %in% HBM_DTXID) # 3944 records

dtxsid_list <- HBM_cyto_bio$dtxsid[1:3945] # 45 compounds
unique(dtxsid_list) # unique DTXSID 

# Save data
save(HBM_cyto_bio, file = '8.4-HBM_Cytotoxicity.RData')

### Clean Cytotoxicity Data
# Load HBM_cytotoxicity data
load('8.4-HBM_Cytotoxicity.RData')
HBM_Cyto <- HBM_cyto_bio #3944 records

# Extract variables from data
HBM_Cyto_1 <- subset(HBM_Cyto, 
                     select = c(ac50, hitcall, dtxsid, maxMeanConc, aeid, spid))

# Add a new variable
library(dplyr)
HBM_Cyto_2 <- HBM_Cyto_1 %>%
  mutate(HIT_CALL = factor(case_when(
    hitcall < 0.9 ~ "Inactive",
    hitcall >= 0.9 ~ "Active"
  )))

## Output Cytotoxicity AC50 value

# Extract AEID, AC50, DTXSID
HBM_Cyto_AC50 <- subset(HBM_Cyto_2,
                        select = c(dtxsid, aeid, spid, ac50)) #3945 records

# 使用 pivot_wider 函数将 aeid 的值变为列名
library(tidyr)
library(dplyr)

# 删除重复值
HBM_Cyto_AC50_unique <- HBM_Cyto_AC50 %>%
  distinct(dtxsid, aeid, ac50, .keep_all = TRUE) #3675 records

# ac50有重复值，因此需要为每个 dtxsid, aeid, spid 组合创建一个唯一的标识符
HBM_Cyto_AC50_unique <- HBM_Cyto_AC50_unique %>%
  group_by(dtxsid, aeid, spid) %>%
  mutate(record_id = row_number()) %>%
  ungroup()

# 使用 pivot_wider 并将 ac50 值展开为单独记录
HBM_Cyto_AC50_Assay <- HBM_Cyto_AC50_unique %>%
  pivot_wider(
    names_from = aeid,
    values_from = ac50,
    names_prefix = "ac50_"
  ) %>%
  arrange(dtxsid, record_id) %>%
  select(-record_id) #376 records, 99 cytotoxicity assays

unique(HBM_Cyto_AC50_Assay$dtxsid) #45 DTXSID

### Select the min_AC50 for each compound across different 'spid'#### 
df <- HBM_Cyto_AC50_Assay %>%
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
      NA  # If the entire column is NA, return NA
    } else {
      min_value <- min(., na.rm = TRUE)
      if (is.infinite(min_value)) NA else min_value
    }
  }))

HBM_Cyto_AC50_Min <- df
save(HBM_Cyto_AC50_Min, file = '8.5-HBM_Cyto_AC50_Min_N=45.RData')

#### Output Cytotoxicity HIT_CALL, calculate %Active assays ####

# Extract aeid, dtxsid, HIT_CALL
HBM_Cyto_HIT <- subset(HBM_Cyto_2,
                       select = c(dtxsid, aeid, spid, HIT_CALL)) #3944 records

# 删除重复值
HBM_Cyto_HIT_unique <- HBM_Cyto_HIT %>%
  distinct(dtxsid, aeid, HIT_CALL, .keep_all = TRUE) #2835 records

# HIT_CALL有重复值，因此需要为每个 dtxsid 和 aeid, spid 组合创建一个唯一的标识符
HBM_Cyto_HIT_unique <- HBM_Cyto_HIT_unique %>%
  group_by(dtxsid, aeid, spid) %>%
  mutate(record_id = row_number()) %>%
  ungroup()

# 使用 pivot_wider 并将 HIT_CALL 值展开为单独记录
HBM_Cyto_HIT_Assay <- HBM_Cyto_HIT_unique %>%
  pivot_wider(
    names_from = aeid,
    values_from = HIT_CALL,
    names_prefix = "HIT_"
  ) %>%
  arrange(dtxsid, record_id) %>%
  select(-record_id) #287 records, 99 cytotoxicity assays

unique(HBM_Cyto_HIT_Assay$dtxsid) #45 DTXSID

### Calculate the percentage of active assays

# 将因子变量转换为数值形式（active 为 1，inactive 为 0）
HBM_Cyto_HIT_Assay_numeric <- HBM_Cyto_HIT_Assay %>%
  mutate(across(starts_with("HIT"), ~ ifelse(. == "Active", 1, 0)))

# 查看数据框结构
str(HBM_Cyto_HIT_Assay_numeric)
# 确认转换后的数据框
print(HBM_Cyto_HIT_Assay_numeric)


# 计算每个 dtxsid 的 active 百分比
HBM_Cyto_HIT_Assay_percentage <- HBM_Cyto_HIT_Assay_numeric %>%
  group_by(dtxsid) %>%
  summarise(
    total_hits = sum(!is.na(across(starts_with("HIT_")))), # 因子变量的总数，忽略NA
    active_count = sum(across(starts_with("HIT_"), ~ sum(. == 1, na.rm = TRUE))), # 计算所有因子变量中的 active 总和
    ratio = (active_count / total_hits) * 100, # 计算百分比
    .groups = 'drop' # 确保正确分组
  )

# 查看计算结果
print(HBM_Cyto_HIT_Assay_percentage) # 45 records

# 合并计算结果, Delete repeated dtxsid record
HBM_Cyto_HIT_Assay_Ratio <- HBM_Cyto_HIT_Assay %>%
  left_join(HBM_Cyto_HIT_Assay_percentage, by = "dtxsid") %>%
  distinct(dtxsid, .keep_all = TRUE) #45 records

# 查看结果
print(HBM_Cyto_HIT_Assay_Ratio)

save(HBM_Cyto_HIT_Assay_Ratio, file = '8.6-HBM_Cyto_%ActiveAssay_N=45.RData')

#### Output Cytoxicity maxMeanConc #####
# Extract AEID, maxMeanConc, DTXSID
HBM_Cyto_Max <- subset(HBM_Cyto_2,
                       select = c(dtxsid, aeid, spid, maxMeanConc)) #3944 records

# 使用 pivot_wider 函数将 aeid 的值变为列名
library(tidyr)
library(dplyr)

# 删除重复值
HBM_Cyto_Max_unique <- HBM_Cyto_Max %>%
  distinct(dtxsid, aeid, maxMeanConc, .keep_all = TRUE) #3709 records

# maxMeanConc有重复值，因此需要为每个 dtxsid 和 aeid, spid 组合创建一个唯一的标识符
HBM_Cyto_Max_unique <- HBM_Cyto_Max_unique %>%
  group_by(dtxsid, aeid, spid) %>%
  mutate(record_id = row_number()) %>%
  ungroup()

# 使用 pivot_wider 并将 maxMeanConc 值展开为单独记录
HBM_Cyto_Max_Assay <- HBM_Cyto_Max_unique %>%
  pivot_wider(
    names_from = aeid,
    values_from = maxMeanConc,
    names_prefix = "max_"
  ) %>%
  arrange(dtxsid, record_id) %>%
  select(-record_id) #338 records, 99 cytotoxicity assays

unique(HBM_Cyto_Max_Assay$dtxsid) #45 DTXSID

library(writexl)
save(HBM_Cyto_Max_Assay, file = '8.7-HBM_Cyto_maxMeanConc_Assay.RData')
