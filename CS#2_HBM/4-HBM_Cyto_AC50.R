
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

#### API-CompTox-ccdR/ctxR package, extract Cytotoxicity assay data ####
library(ctxR)
### API key will be stored as the variable my_key
my_key <- '5e944183-665f-4e1b-9ff0-c7e27b8b09e4'

### Retrieve all assays
assays <- get_all_assays(API_key = my_key) #1499 assays

## Select assays
library(dplyr)

# Select human, liver, cytotoxicity
assays_1 <- assays %>%
  filter(organism == 'human' & 
           intendedTargetFamilySub == 'cytotoxicity') #104 assays

library(writexl)
write_xlsx(assays, '4.1-Human_Cytotoxicity_Assays_Info.xlsx')

# Extract AEID for cytotoxicity and convert to Values
cyto_AEID <- assays_1 %>%
  select(aeid) %>%
  distinct() %>%
  pull(aeid) #104 AEIDs

### Retrieve bioactivity data from cytotoxicity_AEID batch
cyto_bioactivity <- get_bioactivity_details_batch(AEID = cyto_AEID,
                                                  API_key = my_key) #104 cytotoxicity assays

## Merge multiple lists into a database
library(dplyr)
library(purrr)

cyto_bio_merge <- cyto_bioactivity %>%
  map_df(~ .x) # records

save(cyto_bio_merge, file = '4-Human_Cytotoxicity_bio_merge.RData')


##### Extract Cytoxicity data for 76 HBM overlapping compounds ####

## Load Human Cytotoxicity data ##
load('4-Human_Cytotoxicity_bio_merge.RData')

## Load HBM Overlapping compounds N=76 ##
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

dtxsid_list <- HBM_cyto_bio$dtxsid[1:3944] # 45 compounds
unique(dtxsid_list) # unique DTXSID 

# Save data
library(writexl)
write_xlsx(HBM_cyto_bio, '4.2-Overlap_HBM_Cytotoxicity.xlsx')

##### Clean CompTox Cytoxicity Data #####
library(readxl)
HBM_Cyto <- read_excel('4.2-Overlap_HBM_Cytotoxicity.xlsx',
                       sheet = 1) #3944 records

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

#### Output Cytotoxicity AC50 value #####

# Extract AEID, AC50, DTXSID
HBM_Cyto_AC50 <- subset(HBM_Cyto_2,
                          select = c(dtxsid, aeid, spid, ac50)) #3944 records

# 使用 pivot_wider 函数将 aeid 的值变为列名
library(tidyr)
library(dplyr)

# 删除重复值
HBM_Cyto_AC50_unique <- HBM_Cyto_AC50 %>%
  distinct(dtxsid, aeid, ac50, .keep_all = TRUE) #3676 records

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
  select(-record_id) #375 records, 99 cytotoxicity assays

unique(HBM_Cyto_AC50_Assay$dtxsid) #45 DTXSID

library(writexl)
write_xlsx(HBM_Cyto_AC50_Assay, '4.3-Overlap_HBM_Cyto_AC50_Assay.xlsx')


