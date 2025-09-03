
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

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
           tissue == 'liver' &
           intendedTargetFamilySub == 'cytotoxicity') #23 assays

library(writexl)
write_xlsx(assays_1, '3.1-Liver_Cytotoxicity_Assays_Info.xlsx')

# Extract AEID for cytotoxicity and convert to Values
cyto_AEID <- assays_1 %>%
  select(aeid) %>%
  distinct() %>%
  pull(aeid) #23 AEIDs

### Retrieve bioactivity data from cytotoxicity_AEID batch
cyto_bioactivity <- get_bioactivity_details_batch(AEID = cyto_AEID,
                                                  API_key = my_key) #23 cytotoxicity assays

## Merge multiple lists into a database
library(dplyr)
library(purrr)

cyto_bio_merge <- cyto_bioactivity %>%
  map_df(~ .x) #189299 records

save(cyto_bio_merge, file = '3-Human_Liver_Cytotoxicity_bio_merge.RData')

##### Extract Cytoxicity data for 94 liver_related carcinogens ####

## Load liver_related_carcinogens
library(readxl)
liver_carcinogens <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                                sheet = 1)
# Extract DTXID and save as Values
liver_DTXID <- liver_carcinogens %>%
  select(DTXSID) %>%
  distinct() %>%
  pull(DTXSID) #94 DTXSID

## Select liver_related_carcinogens using dtxsid
liver_cyto_bio <- cyto_bio_merge %>%
  filter(dtxsid %in% liver_DTXID) #2649 records

dtxsid_list <- liver_cyto_bio$dtxsid[1:2649]
unique(dtxsid_list) #82 unique DTXSID 

# Save data
library(writexl)
write_xlsx(liver_cyto_bio, '3.2-Overlap_Liver_Cytotoxicity_N=82.xlsx')

##### Clean CompTox Cytoxicity Data #####
library(readxl)
Liver_Cyto <- read_excel('3.2-Overlap_Liver_Cytotoxicity_N=82.xlsx',
                       sheet = 1) #2649 records

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

#### Output Cytotoxicity AC50 value #####

# Extract AEID, AC50, DTXSID
Liver_Cyto_AC50 <- subset(Liver_Cyto_2,
                        select = c(dtxsid, aeid, spid, ac50)) #2649 records

# Use the pivot_wider function to change the aeid value into a column name
library(tidyr)
library(dplyr)

# Remove duplicate values
Liver_Cyto_AC50_unique <- Liver_Cyto_AC50 %>%
  distinct(dtxsid, aeid, ac50, .keep_all = TRUE) #2693 records

# There are duplicate values for ac50, so a unique identifier needs to be created for each dtxsid, aeid, spid combination
Liver_Cyto_AC50_unique <- Liver_Cyto_AC50_unique %>%
  group_by(dtxsid, aeid, spid) %>%
  mutate(record_id = row_number()) %>%
  ungroup()

# Use pivot_wider and expand the ac50 values into separate records
Liver_Cyto_AC50_Assay <- Liver_Cyto_AC50_unique %>%
  pivot_wider(
    names_from = aeid,
    values_from = ac50,
    names_prefix = "ac50_"
  ) %>%
  arrange(dtxsid, record_id) %>%
  select(-record_id) #222 records

unique(Liver_Cyto_AC50_Assay$dtxsid) #82 DTXSID

library(writexl)
write_xlsx(Liver_Cyto_AC50_Assay, '3.3-Overlap_Liver_Cyto_AC50_Assay_N=82.xlsx')




