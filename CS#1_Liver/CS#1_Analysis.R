

####### 1. Overlaping compounds: IARC and NTP_Liver ------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

##### 1.1. Overlap Liver compounds between NTP and IARC list ####
library(readxl)
# Load IARC 2A&2B data
IARC <- read_excel('Chemical_IARC_NTP/IARC_2A&2B.xlsx', sheet = 1)

# Load NTP_Liver data
NTP_Liver <- read_excel('Chemical_IARC_NTP/NTP_liver_gland_kidney_unique.xlsx',
                        sheet = 'liver_unique')

#create list value
IARC_list <- IARC$`CAS No.`[1:369]
liver_list <- NTP_Liver$CASRN[1:184]

#create a list include all list value
venn_list <- list (IARC_2A_2B = IARC_list,
                   NTP_Liver = liver_list) 

#remove NA value from every vector in list
venn_list = purrr::map(venn_list, na.omit)

## Overlap: IARC and NTP_liver
library(VennDiagram)
library(grid)
library(futile.logger)
library(ggvenn)

ggvenn(
  venn_list[c(1,2)], 
  fill_color = c("#0073C2FF", "#EFC000FF"),text_size = 5,
  stroke_size = 0.1, set_name_size = 5, show_percentage = FALSE)

inter <- get.venn.partitions(venn_list[c(1,2)])
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
inter <- subset(inter, select = -..values.. )
inter <- subset(inter, select = -..set.. )

#export intersection results
write.table (inter, '1-Overlap_IARC_NTP_Liver.csv', 
             row.names = FALSE,
             sep = ',',
             quote = FALSE)

#95 overlapping compounds
INC_95_string <- inter$values[1]
INC_95_vector <- unlist(strsplit(INC_95_string, "\\|"))
INC_95_df <- data.frame(CAS_Number = trimws(INC_95_vector))

## Mannually search compounds identifer through CompTox and PubChem according to CAS_Number
# save as Overlap_IARC_NTP.xlsx
#
Liver_95 <- read_excel('Overlap_IARC_NTP_N=95.xlsx',
                       sheet = 'IARC_Liver')

# "Coconut oil diethanolamine condensate" lack INCHIKEY, remove it
Liver_94 <- Liver_95 %>%
  filter(Agent != 'Coconut oil diethanolamine condensate') %>%
  rename(InChIKey = INCHIKEY,
         IARC_Group = Group)

library(writexl)
write_xlsx(Liver_94, '1.1-Overlap_IARC_NTP_N=94.xlsx')

##### 1.2. Extract 94 compounds from Exposome data #####

library(dplyr)
Overlap_Liver <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                            sheet = 1)
## Extract INCHIKEY List
INC_List <- Overlap_Liver$InChIKey[1:94]

# Blood Exposome
load('/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/Blood_Exposome.RData')

# Tissue Exposome
load('/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/Tissue_Exposome.RData')

# Load All_HBM_Merge data
load('/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/All_HBM_Merge_New.RData')

# Load Exposome-Explorer data
load('/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/Exposome_Explorer.RData')

## Extract INCHIKEY List
INC_List <- Overlap_Liver$InChIKey[1:94]

## Extract 94 liver_compounds from HBM and Exposome-Explorer data

# Blood Exposome
Liver_BE <- Blood_Exposome[Blood_Exposome$InChIKey %in% INC_List, ] %>%
  select(InChIKey) %>%
  mutate(Sample = 'Blood',
         Data_Source = 'Blood_Exposome') #78 compounds

# Tissue Exposome
Liver_TE <- Tissue_Exposome[Tissue_Exposome$InChIKey %in% INC_List, ] %>%
  select(InChIKey,Sample_Type) %>%
  rename(Sample = Sample_Type) %>%
  mutate(Data_Source = 'Tissue_Exposome') %>%
  distinct() #78 compounds


# All_HBM
Liver_HBM <- All_HBM_Merge_New[All_HBM_Merge_New$InChIKey %in% INC_List,] %>%
  select(Chemical_Name, InChIKey, CAS_Number, DTXSID, Average_Mass,
         Concentration_Unit, P50, `Geometric Mean`, Mean, P75, P90, P95,
         Sample, Sample_Type) %>%
  # select concentration value
  mutate(Conc = coalesce(P50, `Geometric Mean`, Mean, P75, P90, P95)) %>%
  select(InChIKey, Conc, Concentration_Unit, Sample, Sample_Type) %>%
  rename(Unit = Concentration_Unit)


# Exposome-Explorer
Liver_EE <- Exposome_Explorer[Exposome_Explorer$InChIKey %in% INC_List,] %>%
  select(InChIKey, Geometric.mean, Median, Unit, Biospecimen) %>%
  # select concentration value
  mutate(Conc = coalesce(Median, Geometric.mean)) %>%
  select(InChIKey, Conc, Unit, Biospecimen) %>%
  rename(Sample = Biospecimen) %>%
  mutate(Sample_Type = '')


### Merge all measured data, HBM and EE
Liver_Measured <- bind_rows(Liver_HBM, Liver_EE)

write_xlsx(Liver_Measured, '1.2-Liver_Expo_Conc.xlsx')

### Merge TE, BE, HBM, EE screen results ####
Liver_HBM_1 <- Liver_HBM %>%
  select(InChIKey, Sample) %>%
  mutate(Data_Source = 'HBM') %>%
  mutate(Sample = case_when(
    Sample == "Human (Urine)" ~ "Urine",
    Sample == "Breast Milk" ~ "Breast milk",
    Sample == "hemoglobin adduct" ~ "Hemoglobin adduct",
    Sample == "urine" ~ "Urine",
    Sample == "whole blood" ~ "Blood",
    Sample == "blood" ~ "Blood",
    Sample == "plasma" ~ "Blood plasma",
    Sample == "pooled serum" ~ "Blood serum",
    TRUE ~ Sample  # Keep other unmatched values
  )) %>%
  distinct()

unique(Liver_HBM_1$InChIKey) #27 compounds
unique(Liver_HBM_1$Sample)

Liver_EE_1 <- Liver_EE %>%
  select(InChIKey, Sample) %>%
  mutate(Data_Source = 'Exposome_Explorer') %>%
  mutate(Sample = case_when(
    Sample == "Hemoglobin"  ~ "Blood",
    Sample == "Blood, unspecified" ~ "Blood",
    Sample == "Plasma, unspecified" ~ "Blood plasma",
    Sample == "Serum, fasting" ~ "Blood serum",
    Sample == "Serum, unspecified" ~ "Blood serum",
    Sample == "Whole blood, fasting" ~ "Blood",
    Sample == "Urine, first morning spot" ~ "Urine",
    Sample == "Urine, spot"  ~ "Urine",
    Sample == "Adipose tissue, breast"  ~ "Breast adipose tissue",
    TRUE ~ Sample  # Keep other unmatched values
  ))

unique(Liver_EE_1$Sample)
unique(Liver_EE_1$InChIKey) #8 compounds

# Merge
Liver_merge <- bind_rows(Liver_TE, Liver_BE, Liver_HBM_1, Liver_EE_1)
unique(Liver_merge$InChIKey) #79 compounds

library(dplyr)

Liver_merge_1 <- Liver_merge %>%
  group_by(InChIKey) %>%
  summarise(
    Sample = paste(unique(Sample), collapse = "|"),
    Data_Source = paste(unique(Data_Source), collapse = "|"),
    .groups = 'drop'
  ) #79 compounds

Liver <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                    sheet = 1)

Liver_Screen <- left_join(Liver, Liver_merge_1, by = "InChIKey")

write_xlsx(Liver_Screen, '1.3-Overlap_Liver_Screen_N=94.xlsx')

##### 1.3. Check 94 compounds in Exposome data #####

Liver_BE_inchikey <- Liver_BE %>%
  select(InChIKey)
Liver_TE_inchikey <- Liver_TE %>%
  select(InChIKey)


### Calculate the number of inchikeys that overlap between blood exposome and measured data
Overlap_INC_Count <- 
  length(intersect(Liver_BE$InChIKey, Liver_Measured$InChIKey)) #27 compounds

### Count the number of distinct values in InChIKey
Liver_Measured_INC <- Liver_Measured %>%
  select(InChIKey) %>%
  distinct() #28 compounds

### Count the number of InChIKeys tested
Liver_num_inchikey <- bind_rows(Liver_TE_inchikey,
                                Liver_BE_inchikey,
                                Liver_Measured_INC) %>%
  distinct() ## 79 compounds

# Only in Tissue_Exposome
only_in_TE <- setdiff(Liver_TE$InChIKey, Liver_Measured$InChIKey)
only_in_TE <- setdiff(Liver_TE$InChIKey, Liver_BE$InChIKey)
# 10 compounds in TE are also in BE

# Only in Blood_Exposome
only_in_BE <- setdiff(Liver_BE$InChIKey, Liver_Measured$InChIKey)

# Only in Liver_Measured
only_in_LM <- setdiff(Liver_Measured$InChIKey, Liver_BE$InChIKey)

## Find non-overlap compounds
non_overlap_inchikey <- union(only_in_LM, only_in_BE) #52 compounds

# overlap compounds between BE and Liver_Measured are 79-52=27 compounds

## Find non-overlap compounds in HBM_Measured data
unique_HBM_Measure <- Liver_Measured_INC %>% filter(InChIKey %in% only_in_LM)
# inchikey: QBYJBZPUGVGKQQ-SJJAEHHWNA-N, Aldrin

### Calculate the number of InChIKeys with concentration values
Liver_Conc <- Liver_Measured %>%
  group_by(InChIKey) %>%
  filter(any(!is.na(Conc))) %>%
  select(InChIKey) %>%
  distinct() # 18 compounds

### Count the number of InChIKeys with concentration values of NA
Liver_NA_Conc <- Liver_Measured %>%
  group_by(InChIKey) %>%
  filter(all(is.na(Conc))) %>%
  select(InChIKey) %>%
  distinct() # 10 compounds

### Count the number of InChIKeys with a concentration value of 0
Liver_0_Conc <- Liver_Measured %>%
  group_by(InChIKey) %>%
  filter(all(Conc == 0.0000, na.rm = FALSE)) %>%
  select(InChIKey) %>%
  distinct() # 1 compound

####### 2. Plot_Liver_Conc --------------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

library(readxl)

Liver <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                    sheet = 1)

# Load Liver_Conc data
Liver_Conc <- read_excel('1.2-Liver_Expo_Conc.xlsx',
                         sheet = 1)
unique(Liver_Conc$InChIKey)
##### 2.1. Unit conversion ####

unique(Liver_Conc$Unit)

library(dplyr)

Liver_Conc_1 <- Liver_Conc %>%
  mutate(
    Unit = ifelse(Unit == "ng/mL", "µg/L", Unit),
    Unit = ifelse(Unit == "pmol/g Hb", "pmol/g hemoglobin", Unit),
    Unit = ifelse(Unit == "pmol/g", "pmol/g hemoglobin", Unit),
    
    Unit = ifelse(Unit == "ng/g serum", "ng/g", Unit),
    
    Unit = ifelse(Unit == "ug/g", "µg/g", Unit),
    Unit = ifelse(Unit == "ug/L", "µg/L", Unit),
  ) %>%
  mutate(
    Unit = case_when(
      Unit == "µg/g" & Sample == "Serum, fasting" ~ "µg/g",
      Unit == "µg/g" & Sample == "Breast milk" ~ "µg/g",
      Unit == "µg/g" & Sample == "Adipose tissue, breast" ~ "µg/g",
      Unit == "µg/g" & Sample == "Plasma, unspecified" ~ "µg/g",
      TRUE ~ Unit)  # 
  )

unique(Liver_Conc_1$Unit)

##### 2.2. Sample conversion #####
unique(Liver_Conc$Sample)
unique(Liver_Conc$Sample_Type)

Liver_Conc_2 <- Liver_Conc_1 %>%
  mutate(
    Sample_New = case_when(
      Sample == "Human (Urine)" ~ "Urine",
      Sample == "urine" ~ "Urine",
      Sample == "Urine, spot" ~ "Urine",
      Sample == "Urine, spot" ~ "Urine",
      Sample == "Urine, first morning spot" ~ "Urine",
      
      Sample == "whole blood" ~ "Blood",
      Sample == "blood" ~ "Blood",
      Sample == "Blood, unspecified" ~ "Blood",
      Sample == "Whole blood, fasting" ~ "Blood",
      
      Sample == "plasma"  ~ "Blood plasma",
      Sample == "Plasma, unspecified"  ~ "Blood plasma",
      
      Sample == "Serum" ~ "Blood serum",
      Sample == "pooled serum" ~ "Blood serum",
      Sample == "Serum, fasting" ~ "Blood serum",
      Sample == "Serum, unspecified" ~ "Blood serum",
      
      Sample == "hemoglobin adduct" ~ "Blood",
      Sample == "Hemoglobin" ~ "Blood",
      
      Sample == "Amniotic Fluid" ~ "Amniotic fluid",
      
      Sample == "Adipose tissue, breast" ~ "Breast adipose tissue",
      Sample == "Breast Milk"   ~ "Breast milk",
      
      Sample_Type == "Blood Serum" ~ "Blood serum",
      Sample_Type == "Blood Serum" ~ "Blood serum",
      Sample_Type == "Cord Blood Plasma" ~ "Cord blood plasma",
      
      TRUE ~ Sample)
  )
unique(Liver_Conc_2$Sample_New)

#####

## Merge Sample and Unit variable as Sample_Unit variable
Liver_Conc_3 <- Liver_Conc_2 %>%
  mutate(Sample_Unit = paste(Sample_New, Unit, sep = "_"))

unique(Liver_Conc_3$Sample_Unit)

## Input chemical identifier
Liver_Conc_4 <- left_join(Liver_Conc_3, Liver, by = 'InChIKey') %>%
  rename(Chemical_Name = PREFERRED_NAME)

unique(Liver_Conc_4$InChIKey) #28 compounds

write_xlsx(Liver_Conc, '2.1-Liver_Conc_Unit_N=28.xlsx')

## Calculate P50, P5, P95 for each compound
Liver_Conc_5 <- Liver_Conc_4 %>%
  filter(!is.na(Conc)) %>%
  group_by(InChIKey, Chemical_Name, Sample_Unit) %>%
  summarize(
    P5 = quantile(Conc, 0.05, na.rm = TRUE),
    P50 = quantile(Conc, 0.50, na.rm = TRUE),
    P95 = quantile(Conc, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # delete conc = 0.0000 records
  filter(P50 != 0)

unique(Liver_Conc_5$Sample_Unit)

write_xlsx(Liver_Conc_5, '2.2-Liver_Conc_Unit_P5.P50.P95.xlsx')

#### Creating a graph: Adding group labels ####

library(tidyr)
library(dplyr)

# Make sure all compounds are visible on the y-axis
Liver_Conc_5 <- Liver_Conc_5 %>%
  mutate(Chemical_Name = factor(Chemical_Name, levels = unique(Chemical_Name))) 


library(ggplot2)
library(scales)

custom_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
  "#0000FF", "#A65628", "#F781BF", "#808080", "#045F5F", 
  "#800517", "#8DA0CB", "#EDDA74", "#A6D854", "#ADD8E6",
  "#FF00FF", "#FFC0CB", "#7FFFD4", "#00FFFF")

ggplot(Liver_Conc_5, aes(y = reorder(Chemical_Name, P5))) +
  # Scatter plot with custom shape based on Sample_Unit
  geom_point(aes(x = P50, 
                 color = Sample_Unit, 
                 shape = Sample_Unit),  # Map shape to Sample_Unit
             size = 3,  alpha = 0.7, , stroke = 0.7) +
  # Error bars (no shape aesthetic here)
  geom_errorbar(aes(xmin = P5, xmax = P95,
                    colour = Sample_Unit),
                width = 0.3, linewidth = 0.8) +
  # Set x-axis to logarithmic scale
  scale_x_log10(limits = c(1e-4, 1e3),
                breaks = 10^seq(-4, 3, by = 1),
                labels = trans_format("log10", math_format(10^.x))) +
  
  # Custom shape mapping, blood and serum are set to 1, urine is set to 2
  scale_shape_manual(values = c("Blood_µg/L" = 1,
                                "Blood serum_ng/g" = 1,
                                "Blood serum_ng/g lipid" = 1, 
                                "Blood serum_µg/L" = 1,
                                "Blood serum_µg/g lipid" = 1,
                                "Blood_pmol/g hemoglobin" = 1,
                                "Blood_ng/L" = 1,
                                "Blood_pg/mL" = 1,
                                "Blood plasma_ng/g lipid" = 1,
                                "Blood plasma_µg/L"  = 1,
                                
                                "Breast milk_ng/g" = 4,
                                "Breast milk_µg/L" = 4,
                                
                                "Urine_µg/L" = 2, 
                                "Urine_µg/g" =2,
                                "Urine_µg/g creatinine" = 2,
                                
                                "Cord blood plasma_µg/L" = 3,
                                "Cord blood plasma_µg/g lipid" = 3,
                                
                                "Breast adipose tissue_ng/g" = 7,
                                "Semen_µg/L" = 9
  )) +
  # set color
  scale_color_manual(values = custom_colors) +
  
  # Add labels and theme
  labs(
    x = "Concentration (log scale)",
    y = ''
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11),  # Adjust x-axis label size
    axis.text.y = element_text(size = 10),  # Adjust y-axis label size
    legend.position = "right"       # Place legend on the right 
  )



####### 3. Liver_Cyto_AC50 -----------------------------------------------------

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


####### 4.Liver_HTTK_HExp_SEEM3-------------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

library(readxl)

Che_list <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                       sheet = 1)

### HTTK model to calculate Css ###
library(httk)
library(readxl)
library(ggplot2)
library(dplyr)

### Step 1: Loading In Vitro Screening Data ####

load('invitrodb_3_5_mc5.Rdata')

# Built new subset of chemicals with human HTTK data
toxcast.httk <- subset(mc5, dsstox_substance_id %in% get_cheminfo(
  info='DTXSID', suppress.messages = TRUE))

# Extract 'DTXSID' from our chemcial list
my.chems <- Che_list$DTXSID #76

# Extract cooresponding chemicals from Toxcast
example.toxcast <- as.data.frame(mc5[mc5$dsstox_substance_id %in% my.chems,]) #51074 records

# pick 9 variables
example.toxcast <- example.toxcast[, c("chnm",
                                       "dsstox_substance_id",
                                       "spid",
                                       "hitc",
                                       "modl",
                                       "aeid",
                                       "modl_ga",
                                       "modl_ac10",
                                       "modl_acc")]

# reduce precision to decrease space:
cols <- c("modl_ga", "modl_ac10", "modl_acc")
for (this.col in cols)
  example.toxcast[, this.col] = signif(example.toxcast[, this.col], 3)

#Each chemical might have multiple assays for which there is a hit. 
#Each chemical may have been run multiple times (experimental replicates), 
#so that there may be multiple samples from the same chemical for the same assay. 
#Different samples are indicated by the column “spid” (sample id).

knitr::kable(head(example.toxcast), caption = "ToxCast In Vitro Bioactivity Data",
             floating.environment="sidewaystable")

# Find the in vitro concentrations for compounds
toxcast.table <- NULL 
for (this.id in unique(example.toxcast$dsstox_substance_id))
{
  this.subset <- subset(example.toxcast, dsstox_substance_id == this.id)
  these.hits <- subset(this.subset, hitc==1)
  if (dim(these.hits)[1]>0){
    this.row <- data.frame(Compound=as.character(this.subset[1,"chnm"]),
                           DTXSID=this.id,
                           Total.Assays = dim(this.subset)[1],
                           Unique.Assays = length(unique(this.subset$aeid)),
                           Total.Hits = dim(these.hits)[1],
                           Unique.Hits = length(unique(these.hits$aeid)),
                           Low.AC50 = signif(min(these.hits$modl_ga),3),
                           Low.AC10 = signif(min(these.hits$modl_ac10),3),
                           Low.ACC = signif(min(these.hits$modl_acc),3),
                           Q10.AC50 = signif(quantile(these.hits$modl_ga,probs=0.1),3),
                           Q10.AC10 = signif(quantile(these.hits$modl_ac10,probs=0.1),3),
                           Q10.ACC = signif(quantile(these.hits$modl_acc,probs=0.1),3),
                           Med.AC50 = signif(median(these.hits$modl_ga),3),
                           Med.AC10 = signif(median(these.hits$modl_ac10),3),
                           Med.ACC = signif(median(these.hits$modl_acc),3),
                           Q90.AC50 = signif(quantile(these.hits$modl_ga,probs=0.9),3),
                           Q90.AC10 = signif(quantile(these.hits$modl_ac10,probs=0.9),3),
                           Q90.ACC = signif(quantile(these.hits$modl_acc,probs=0.9),3)
    )
    toxcast.table <- rbind(toxcast.table, this.row)
  }
}
rownames(toxcast.table) <- seq(1,dim(toxcast.table)[1]) #39
knitr::kable(head(toxcast.table[,1:6]), caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")


# Save HBM_HTTK_toxcast.table
library(writexl)
write_xlsx(toxcast.table, '4.1-Liver_HTTK_toxcast_table.xlsx')

### Step 2: HTTK deduce Css_Plasma with in vitro ####

#Steady State Plasma Concentration Css 
#calc_mc_css(): predict Css resulting from an ongoing 1 mg/kg/day exposure.

# 
for (this.id in unique(toxcast.table$DTXSID)) {
  # 
  if (this.id %in% get_cheminfo(info = "dtxsid", suppress.messages = TRUE)) {
    # 
    set.seed(12345)
    
    # 
    css_results <- calc_mc_css(dtxsid = this.id,
                               output.units = "uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # httk default sampling number
                               suppress.messages = TRUE)
    
    # 
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_5th"]  <- css_results[3]
    
    # 
    toxcast.table[toxcast.table$DTXSID == this.id, "Css.Type"] <- "in vitro"
  }
}

#Let’s look at just the relevant columns from our table. 
#Any Css values of “NA” (not a number) indicate that we don’t currently have in vitro TK data for those chemicals:

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th', "Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")


### Step 3: HTTK deduce Css_Plasma with QSPR Predictions ####

#For chemicals where no in vitro toxicokinetic data are available, 
#we can sometimes load predictions from in silico models. 

load_sipes2017()

library(httk)

# 
if (!"Css_95th" %in% colnames(toxcast.table)) {
  toxcast.table$Css_95th <- NA
}
if (!"Css_50th" %in% colnames(toxcast.table)) {
  toxcast.table$Css_50th <- NA
}
if (!"Css_5th" %in% colnames(toxcast.table)) {
  toxcast.table$Css_5th <- NA
}
if (!"Css.Type" %in% colnames(toxcast.table)) {
  toxcast.table$Css.Type <- NA
}

# 
for (this.id in unique(toxcast.table$DTXSID)) {
  
  # 
  if (this.id %in% get_cheminfo(info="dtxsid", suppress.messages=TRUE) &
      is.na(toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"])) {
    
    # 
    set.seed(12345)
    
    # 
    css_results <- calc_mc_css(dtxsid=this.id,
                               output.units="uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  
                               suppress.messages=TRUE)
    
    # 
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_5th"] <- css_results[3]
    
    # 
    toxcast.table[toxcast.table$DTXSID==this.id, "Css.Type"] <- "in silico"
  }
}

#Take another look at our table – many of the gaps are now filled in.

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th',"Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")

sum(!is.na(toxcast.table$Css_50th))

## Output
Liver_Css_Plsama = toxcast.table
Liver_Css_Plsama <- Liver_Css_Plsama %>%
  select(DTXSID, Css_95th, Css_50th, Css_5th, Css.Type) %>%
  dplyr::rename(Css_plasma_95th = Css_95th,
                Css_plasma_50th = Css_50th,
                Css_plasma_5th = Css_5th,
                Css.Type_Plasma = Css.Type)

### Deduce Blood Css ####

### Step 2.1: HTTK deduce Css_Blood with in vitro ####

#Steady State Plasma Concentration Css 
#calc_mc_css(): predict Css resulting from an ongoing 1 mg/kg/day exposure.

# 
for (this.id in unique(toxcast.table$DTXSID)) {
  # 
  if (this.id %in% get_cheminfo(info = "dtxsid", suppress.messages = TRUE)) {
    # 
    set.seed(12345)
    
    # 
    css_results <- calc_mc_css(dtxsid = this.id,
                               output.units = "uM",
                               concentration = 'blood',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # httk 默认采样数
                               suppress.messages = TRUE)
    
    # 
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_5th"]  <- css_results[3]
    
    # 
    toxcast.table[toxcast.table$DTXSID == this.id, "Css.Type"] <- "in vitro"
  }
}

#Let’s look at just the relevant columns from our table. 
#Any Css values of “NA” (not a number) indicate that we don’t currently have in vitro TK data for those chemicals:

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th', "Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")


### Step 3.1: HTTK deduce Css_Blood with QSPR Predictions ####

#For chemicals where no in vitro toxicokinetic data are available, 
#we can sometimes load predictions from in silico models. 

load_sipes2017()

library(httk)

# 
if (!"Css_95th" %in% colnames(toxcast.table)) {
  toxcast.table$Css_95th <- NA
}
if (!"Css_50th" %in% colnames(toxcast.table)) {
  toxcast.table$Css_50th <- NA
}
if (!"Css_5th" %in% colnames(toxcast.table)) {
  toxcast.table$Css_5th <- NA
}
if (!"Css.Type" %in% colnames(toxcast.table)) {
  toxcast.table$Css.Type <- NA
}

# 
for (this.id in unique(toxcast.table$DTXSID)) {
  
  # 
  if (this.id %in% get_cheminfo(info="dtxsid", suppress.messages=TRUE) &
      is.na(toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"])) {
    
    # 
    set.seed(12345)
    
    # 
    css_results <- calc_mc_css(dtxsid=this.id,
                               output.units="uM",
                               concentration = 'blood',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  
                               suppress.messages=TRUE)
    
    # 
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_5th"] <- css_results[3]
    
    # 
    toxcast.table[toxcast.table$DTXSID==this.id, "Css.Type"] <- "in silico"
  }
}

#Take another look at our table – many of the gaps are now filled in.

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th',"Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")

sum(!is.na(toxcast.table$Css_50th))

## Output
Liver_Css_Blood = toxcast.table
Liver_Css_Blood <- Liver_Css_Blood %>%
  select(DTXSID, Css_95th, Css_50th, Css_5th, Css.Type) %>%
  dplyr::rename(Css_blood_95th = Css_95th,
                Css_blood_50th = Css_50th,
                Css_blood_5th = Css_5th,
                Css.Type_Blood = Css.Type)

##### Merge plasma and blood Css ####
Liver_Css <- left_join(Liver_Css_Blood, Liver_Css_Plsama, by = 'DTXSID')

### Extract predicted internal blood from HExpPredict ####
load('/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/Internal_Prediction.RData')
Predict_Blood <- Internal_Prediction %>%
  select(DTXSID, P5, P50, P95) %>%
  dplyr::rename(CB_P5_uM = P5,
                CB_P50_uM = P50,
                CB_P95_uM = P95) 

Liver_Predict_Blood <- left_join(Che_list, Predict_Blood, by = 'DTXSID')

### Extract predicted external exposure from SEEM3 ####
load('/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/External_Prediction.RData')
Predict_External <- External_Prediction %>%
  select(DTXSID, P5, P50, P95, Concentration_Unit) %>%
  dplyr::rename(Ex_P5 = P5,
                Ex_P50 = P50,
                Ex_P95 = P95,
                Ex_Unit = Concentration_Unit)

Liver_Predict <- left_join(Liver_Predict_Blood, Predict_External, 
                           by = 'DTXSID')

##### Merge ####
Pred_In_Ex <- left_join(Liver_Predict, Liver_Css_Blood, by = 'DTXSID')
Pred_In_Ex_1 <- left_join(Pred_In_Ex, Liver_Css_Plsama, by = 'DTXSID')

library(writexl)
write_xlsx(Pred_In_Ex_1, '4.2-Overlap_Liver_Prediction.xlsx')

####### 5.Plot_Liver_Pred_Expo_AC50---------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

## Load Liver_Cytotoxicity AC50 data ####

Cyto <- read_excel('3.3-Overlap_Liver_Cyto_AC50_Assay_N=82.xlsx',
                   sheet = 1)

library(tidyr)

Cyto_1 <- Cyto %>%
  select(-spid, -ac50_1971, -ac50_2031, -ac50_2450) %>%
  rename(DTXSID = dtxsid) %>%
  pivot_longer(cols = starts_with("ac50_"), 
               names_to = "ac50_type", 
               values_to = "ac50_value",
               values_drop_na = TRUE) 

unique(Cyto_1$DTXSID) #82 compounds

## Load Liver_Prediction Data ####
Liver_Prediction <- read_excel('4.2-Overlap_Liver_Prediction.xlsx',
                               sheet = 1) %>%
  select(DTXSID, Css_blood_5th, Css_blood_50th, Css_blood_95th, Css.Type_Blood,
         CB_P5_uM, CB_P50_uM, CB_P95_uM,
         Ex_P5, Ex_P50, Ex_P95)

## Load Liver_Expo Data ####
Liver_Conc <- read_excel('2.2-Liver_Conc_Unit_P5.P50.P95.xlsx',
                         sheet = 1) %>%
  mutate(Unit = sub(".*_", "", Sample_Unit),
         Sample = sub("_[^_]+$", "", Sample_Unit)) %>%
  select(-Chemical_Name)

# Input compound AVERAGE_MASS
Liver <- read_excel('Overlap_IARC_NTP_N=95.xlsx',
                    sheet = "IARC_Liver") %>%
  select(INCHIKEY, AVERAGE_MASS) %>%
  rename(InChIKey = INCHIKEY)

Che_List <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                       sheet = 1) %>%
  select(InChIKey, DTXSID, PREFERRED_NAME) %>%
  rename(Chemical_Name = PREFERRED_NAME)

Che_List <- left_join(Che_List, Liver, by = 'InChIKey')

Liver_Conc_1 <- left_join(Che_List, Liver_Conc, by = 'InChIKey')


# Create Sample_Detail
unique(Liver_Conc$Unit)
unique(Liver_Conc$Sample)

# Convert unit "µg/L", "ng/L", "pg/mL" to "µmol/L"
Liver_Conc_2 <- Liver_Conc_1 %>%
  mutate(
    Conc_P5_uM = case_when(
      Unit == "µg/L"  ~ P5 / AVERAGE_MASS,
      Unit == "ng/L"  ~ (P5 * 0.001) / AVERAGE_MASS,
      Unit == "pg/mL" ~ (P5 * 0.001) / AVERAGE_MASS,
      TRUE ~ NA_real_
    ),
    Conc_P50_uM = case_when(
      Unit == "µg/L"  ~ P50 / AVERAGE_MASS,
      Unit == "ng/L"  ~ (P50 * 0.001) / AVERAGE_MASS,
      Unit == "pg/mL" ~ (P50 * 0.001) / AVERAGE_MASS,
      TRUE ~ NA_real_
    ),
    Conc_P95_uM = case_when(
      Unit == "µg/L"  ~ P95 / AVERAGE_MASS,
      Unit == "ng/L"  ~ (P95 * 0.001) / AVERAGE_MASS,
      Unit == "pg/mL" ~ (P95 * 0.001) / AVERAGE_MASS,
      TRUE ~ NA_real_
    ),
    Unit_uM = case_when(
      Unit %in% c("µg/L", "ng/L", "pg/mL") ~ "µmol/L",
      TRUE ~ NA_character_
    )
  )

## Save
library(writexl)
write_xlsx(Liver_Conc_2, '5-Overlap_Liver_Conc_uM_Unit_P5.P50.P95.xlsx')

##### Merge HBM_Conc, HBM_HBM_Prediction, Cyto_AC50 ####

df <- left_join(Liver_Prediction, Liver_Conc_2,  by = 'DTXSID')

df_1 <- left_join(df, Cyto_1, by = 'DTXSID')

unique(df_1$DTXSID) #94 compounds

Liver_Pred_Expo_CytoAC50 <- df_1
save(Liver_Pred_Expo_CytoAC50, file = '5-Liver_Pred_Expo_CytoAC50_N=94.RData')


#### Plot: Cyto_AC50, CB, Css, Expo, SEEM3 #####

load('5-Liver_Pred_Expo_CytoAC50_N=94.RData')
df <- Liver_Pred_Expo_CytoAC50

unique(df$Sample)

library(ggplot2)
library(patchwork)

## Plot 1: CB, Css, Exposome (Conc_uM) ####
plot1 <- ggplot(df, aes(y = reorder(Chemical_Name, Ex_P50))) +  
  
  # Add ac50 points
  geom_point(aes(x = ac50_value, color = "Cytotoxicity_AC50"), 
             size = 1, shape = 3, alpha = 0.7) +
  
  # Add CB error bars to show the range from 5th to 95th percentiles
  geom_errorbar(aes(xmin = CB_P5_uM, 
                    xmax = CB_P95_uM, color = "CB_95%CI"), 
                width = 0.3) +
  # Add CB median point
  geom_point(aes(x = CB_P50_uM, color = "CB_Median"), 
             size = 1.5, shape = 16) +
  
  # Added Css_blood error bars to represent the 5th to 95th percentile range
  geom_errorbar(aes(xmin = Css_blood_5th, 
                    xmax = Css_blood_95th, color = "Css_blood_95%CI"), 
                width = 0.3) +
  
  # Add Css_blood median point
  geom_point(aes(x = Css_blood_50th, color = "Css_blood_Median"), 
             size = 1.5, shape = 16, alpha = 0.7) +
  
  # Add Conc Error Bars
  geom_errorbar(aes(xmin = Conc_P5_uM, 
                    xmax = Conc_P95_uM, color = "Measured_Conc_95%CI"), 
                width = 0.3) +
  
  # Add Conc median point and group different shapes according to Conc_Sample
  geom_point(aes(x = Conc_P50_uM, color = "Measured_Conc_Median", 
                 shape = Sample), size = 1.5, alpha = 0.7) +
  
  # Set the x-axis to logarithmic scale
  scale_x_log10(limits = c(1e-6, 1e5),
                breaks = 10^seq(-6, 5, by = 1),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  
  # Add tags
  labs(x = "Concentration (µM)", color = "Legend", shape = "Sample Type") +
  
  # Custom Colors
  scale_color_manual(values = c("Cytotoxicity_AC50" = '#FFA500',
                                "Css_blood_Median" = "red", 
                                'Css_blood_95%CI'= 'red',
                                "CB_Median" = "blue", "CB_95%CI" = "blue",
                                "Measured_Conc_Median" = "#4E9258", 
                                'Measured_Conc_95%CI'= '#4E9258')) +
  
  # Custom shapes
  scale_shape_manual(values = c("Urine" = 17,  
                                "Blood" = 16,  
                                "Blood plasma" = 1,  
                                "Blood serum" = 5,
                                "Cord blood plasma" = 6,
                                "Plasma" = 1, 
                                "Breast milk" = 4,
                                "Breast adipose tissue" = 9
  )) +  
  
  # Remove the y-axis label and ticks
  theme_minimal() +
  theme(axis.title.y = element_blank(),  
        axis.text.y = element_blank(),   
        axis.ticks.y = element_blank(),  
        axis.text.x = element_text(size = 7))

print(plot1)


## Plot 2: External Prediction SEEM3 ####
plot2 <- ggplot(df, aes(y = reorder(Chemical_Name, Ex_P50))) +  
  
  # Added SEEM3 error bars to represent the range of 5th to 95th percentiles
  geom_errorbar(aes(xmin = Ex_P5, 
                    xmax = Ex_P95, color = "Exposure_95%CI"), 
                width = 0.3) +
  
  # Add SEEM3 median point
  geom_point(aes(x = Ex_P50, color = "Exposure_Median"), 
             size = 1.5, shape = 0) +  # Use different colors and shapes
  
  # Set the x-axis to logarithmic scale
  scale_x_log10(limits = c(1e-17, 1e1),
                breaks = 10^seq(-17, 1, by = 2),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  
  # Add tags
  labs(x = "Concentration (mg/kg-bw/day)", color = "Legend") +
  
  # Custom Colors
  scale_color_manual(values = c("Exposure_Median" = "purple", 
                                'Exposure_95%CI'= 'purple')) +
  
  # Keep y-axis label and ticks
  theme_minimal() +
  theme(axis.title.y = element_blank(),  
        axis.text.y = element_text(size = 7.25),  
        axis.ticks.y = element_line(),      
        axis.text.x = element_text(size = 7.25))

print(plot2)


## Merge two charts, align them horizontally, and share a common y-axis ####
combined_plot <- plot2 + plot1 + 
  plot_layout(ncol = 2, guides = "collect") & theme(legend.position = 'right')

# 
combined_plot

#### 5.1 Overlap_Liver_Sample ####

library(readxl)
library(dplyr)
df <- read_excel('5-Overlap_Liver_Conc_uM_Unit_P5.P50.P95.xlsx',
                 sheet = 1)
unique(df$InChIKey) #94 compounds

df_1 <- df %>%
  mutate(
    Detection = case_when(
      !is.na(Sample) ~ "Yes",
      TRUE ~ "No"
    ),
    Available_Conc_uM = case_when(
      !is.na(Sample) & !is.na(Conc_P50_uM) ~ "Yes",
      TRUE ~ "No"
    )
  )

df_2 <- df_1 %>%
  select(InChIKey, Sample, Detection, Available_Conc_uM) %>%
  mutate(
    Available_Conc_uM_1 = case_when(
      Detection == "Yes" ~ paste0(Available_Conc_uM, ",", Sample),
      TRUE ~ Available_Conc_uM
    )
  )
unique(df_2$InChIKey)

library(stringr)

df_2 <- df_1 %>%
  select(InChIKey, Sample, Detection, Available_Conc_uM) %>%
  group_by(InChIKey) %>%
  mutate(
    Sample_1 = if(any(Detection == "Yes")) {
      str_c(unique(Sample), collapse = ",")
    } else {
      NA_character_  
    }
  ) %>%
  ungroup()

df_3 <- df_2 %>%
  group_by(InChIKey) %>%
  mutate(
    Available_Conc_uM_1 = case_when(
      Available_Conc_uM == "Yes" ~ {
        
        samples_concat <- str_c(unique(Sample[Available_Conc_uM == "Yes"]), collapse = ",")
        
        str_c("Yes", samples_concat, sep = ",")
      },
      Available_Conc_uM == "No" ~ NA_character_,
      TRUE ~ NA_character_  
    )
  ) %>%
  ungroup()

df_4 <- df_3 %>%
  select(-Sample, -Available_Conc_uM) %>%
  distinct()

unique(df_4$InChIKey) #94 compounds


df_5 <- df_4 %>%
  group_by(InChIKey, Sample_1) %>%    
  filter(
    
    !(Detection == "Yes" & n() > 1 & is.na(Available_Conc_uM_1))
  ) %>%
  ungroup() %>%
  rename(Sample = Sample_1, Available_Conc_uM = Available_Conc_uM_1)

library(writexl)
write_xlsx(df_5, '5.1-Overlap_Liver_Sample.xlsx')


####### 6.

####### 6.Liver_LD50_ChemPro_CTV_Cytotoxicity-----------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

#### 6.1. ChemPro and LD50 data for HBM 76 compounds ####
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

##### 6.2. Conditional Toxicity Value (CTV) #####
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

##### 6.3. CompTox Cytotoxicity Data #####

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






####### 7.ToxPi_Liver_Data_Preparation-----------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

library(readxl)
library(dplyr)
##### 7.1. Classify IARC evidence and NTP evidence level #####
## IARC
Liver <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                    sheet = 1)

# 2A=2, 2B=1
Liver_IARC <- Liver %>%
  mutate(IARC_Group = case_when(
    IARC_Group == '2A' ~ 2,
    IARC_Group == '2B' ~ 1))

## NTP evidence
# Positive=3, Clear evidence=2, Some evidence=1
All_NTP_Liver <- read_excel('Chemical_IARC_NTP/NTP_raw_data/2024-06-26-site_data_liver.xlsx',
                            sheet = 'Worksheet') #386 records

All_NTP_Liver_1 <- All_NTP_Liver %>%
  select(CASRN, `Level of Evidence`) %>%
  distinct(CASRN, `Level of Evidence`) %>%
  mutate( `Level of Evidence` = case_when(
    `Level of Evidence` == 'Positive' ~ 3,
    `Level of Evidence` == 'Clear Evidence' ~ 2,
    `Level of Evidence` == 'Some Evidence' ~ 1)) %>%
  group_by(CASRN) %>%
  summarize(across(`Level of Evidence`,  ~ mean(.x, na.rm = TRUE))) %>% #For one compound has multiple evidence level, calculating their mean value
  rename(NTP_Group = `Level of Evidence`)

Liver_NTP_Evi <- All_NTP_Liver_1 %>%
  filter(CASRN %in% Liver$CASRN) #94 compounds

## Merge IARC and NTP evidence
Liver_IARC_NTP <- left_join(Liver_IARC, Liver_NTP_Evi,
                            by = 'CASRN')

write_xlsx(Liver_IARC_NTP, '7.1-Liver_IARC_NTP_N=94.xlsx')

Liver_IARC_NTP_1 <- Liver_IARC_NTP %>%
  select(DTXSID, IARC_Group, NTP_Group)

### 7.2. LD50 Prediction data ####
load('6.1-Liver_LD50_Pred_N=79.RData') #94 compounds
LD50_Pred <- LD50_Pred %>%
  rename(LD50 = `ORAL_RAT_LD50_MOL/KG_TEST_PRED`)
unique(LD50_Pred$DTXSID)

### 7.3. Chemical Properties data ####
load('6.2-Liver_Kow_Koa_Pred_N=88.RData')

# Filter data, OPERA 1st option, EPISUITE 2nd option
# 
Kow_Koa <- Kow_Koa %>%
  filter(SOURCE == 'OPERA 2.6') %>%
  select(DTXSID, NAME, VALUE)

library(tidyr)

Kow_Koa_1 <- Kow_Koa %>%
  
  mutate(VALUE = as.numeric(VALUE)) %>%
  pivot_wider(
    names_from = NAME,  
    values_from = VALUE  
  ) %>% 
  rename(LogKoa = `LogKoa: Octanol-Air`,
         LogKow = `LogKow: Octanol-Water`)#87 compounds

unique(Kow_Koa_1$DTXSID)

### 7.4. CTV Pre+Exp data ####
load('6.3-Liver_CTV_Merge.RData')

# 
CTV_1 <- CTV %>%
  select(-Unit, -Type) %>%
  pivot_wider(names_from = Tox_parameter, values_from = Value) %>%
  select(-PREFERRED_NAME) #63 compounds

CTV_1 <- CTV_1 %>%
  mutate(across(2:9, ~ as.numeric(.)))

### 7.5. Cyto_AC50 ####
load('6.5-Liver_Cyto_AC50_Min_N=82.RData')

AC50_Min <- Liver_Cyto_AC50_Min %>%
  select(-ac50_1971, -ac50_2031, -ac50_2450) %>%
  rename(DTXSID = dtxsid) #82 compounds

### 7.6. Cyto_%Activity Assay ####
load('6.6-Liver_Cyto_%ActiveAssay_N=82.RData')

HIT_Assay_Ratio <- Liver_Cyto_HIT_Assay_Ratio %>%
  select(dtxsid, ratio) %>%
  rename(HIT_Ratio = ratio,
         DTXSID = dtxsid) #82 compounds

##### Merge all data ####

Liver_Tox <- Liver %>%
  left_join(Liver_IARC_NTP_1, by = 'DTXSID') %>%
  left_join(Kow_Koa_1, by = 'DTXSID') %>%
  left_join(LD50_Pred, by = 'DTXSID') %>%
  left_join(CTV_1, by = 'CASRN') %>%
  left_join(HIT_Assay_Ratio, by = 'DTXSID') %>%
  left_join(AC50_Min, by = 'DTXSID')

save(Liver_Tox, file = '7-Liver_Tox_N=94.RData')


####### 8.ToxPi_Plot_Liver------------------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

#### Toxpi_score: Liver
load('7-Liver_Tox_N=94.RData')

###### Toxpi model conctruction ####
library(toxpiR)

# View all functions
lsf.str("package:toxpiR")

## Create slice with transformation functions
# Slice 1: LogKow, LogKoa, log10(x)-log10(min)
# Slice 2: Cyto_AC50, # -log10(x)+log10(max)
# Slice 3: Cyto_%Active, # linear 
# Slice 4: RfD, # -log10(x)+log10(max) 
# Slice 5: BMD, # -log10(x)+log10(max)
# Slice 6: BMDL, # -log10(x)+log10(max)
# Slice 7: NOAEL, # -log10(x)+log10(max)
# Slice 8: CPV, # log10(x)-log10(min)
# Slice 9: OSF, # log10(x)-log10(min)
# Slice 10: IUR, # log10(x)-log10(min) 
# Slice 11: RfC, # -log10(x)+log10(max) 
# Slice 12: LD50, # -log10(x)+log10(max)
# Slice 13: NTP_Group (In vivo evidence), # linear
# Slice 14: IARC, #linear

# 
ac50_columns <- grep("^ac50_", names(Liver_Tox), value = TRUE)

f.slices <- TxpSliceList(
  Slice1 =  TxpSlice(c('LogKoa','LogKow')),
  Slice2 = TxpSlice(ac50_columns),
  Slice3 = TxpSlice('HIT_Ratio'),
  Slice4 = TxpSlice('RfD'),
  Slice5 = TxpSlice('BMD'),
  Slice6 = TxpSlice('BMDL'),
  Slice7 = TxpSlice('NOAEL'),
  Slice8 = TxpSlice('CPV'),
  Slice9 = TxpSlice('OSF'),
  Slice10 = TxpSlice('IUR'),
  Slice11 = TxpSlice('RfC'),
  Slice12 = TxpSlice('LD50'),
  Slice13 = TxpSlice('NTP_Group'),
  Slice14 = TxpSlice('IARC_Group.y')
)

## Goal - Create ToxPi model.
# All Slices, weight = 1.

## Method#2: Mean imputation instead of missing value
# final.trance have missing value, select mean imputation
# Define the mean imputation function
mean_imputation <- function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}
linear <- function(x) {
  # Example linear transformation
  return(x)
}

# Define the transformation functions list

final.trans <- TxpTransFuncList(
  f1 = function(x) { 
    x <- mean_imputation(x)
    log10(x) - log10(min(x))
  },
  f2 = function(x) {
    # Perform mean imputation
    x <- mean_imputation(x)
    # Apply the transformation
    -log10(x) + log10(max(x))
  }, 
  f3 = linear,              
  f4 = function(x) { 
    x <- mean_imputation(x) 
    -log10(x) + log10(max(x))
  }, 
  f5 = function(x) { 
    x <- mean_imputation(x)
    -log10(x) + log10(max(x))
  }, 
  f6 = function(x) { 
    x <- mean_imputation(x)
    -log10(x) + log10(max(x))
  }, 
  f7 = function(x) { 
    x <- mean_imputation(x)
    -log10(x) + log10(max(x))
  }, 
  f8 = function(x) { 
    x <- mean_imputation(x)
    log10(x) - log10(min(x))
  }, 
  f9 = function(x) { 
    x <- mean_imputation(x)
    log10(x) - log10(min(x))
  }, 
  f10 = function(x) { 
    x <- mean_imputation(x)
    log10(x) - log10(min(x))
  },
  f11 = function(x) { 
    x <- mean_imputation(x) 
    -log10(x) + log10(max(x))
  }, 
  f12 = function(x) { 
    x <- mean_imputation(x)
    -log10(x) + log10(max(x))
  },
  f13 = linear,
  f14 = linear
)


f.model <- TxpModel(txpSlices = f.slices,
                    txpWeights = c (1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                    txpTransFuncs = final.trans)

## Calculate ToxPi scores

f.results <- txpCalculateScores(model = f.model,
                                input = Liver_Tox,
                                id.var = 'PREFERRED_NAME.x')
# ToxPi scores
toxpi_scores <- txpSliceScores(f.results)
print(toxpi_scores)

#### Pie Plotting ####

## Setting colors
library(RColorBrewer)
# Get the 12 colors from the "Paired" palette
colors <- brewer.pal(12, "Paired")
# Add two additional distinct colors
colors <- c(colors, "#FF69B4", "#8A2BE2")

## Draw plot
library(grid)
plot(f.results,
     showScore = TRUE,
     fills = colors # fill colors
)

#### Edit Pie Plotting ####
# 
grid.ls()
# 
all_grobs <- grid.ls(print = FALSE)


## 
label_grobs <- grep("label", all_grobs$name, value = TRUE)

for (label in label_grobs) {
  grid.edit(label, gp = gpar(fontsize = 6))
}

## Edit 'Toxpi Scores' characteristics
for (s in sprintf("pie-%d-radSum", 1:94)) {
  grid.edit(s, gp = gpar(fontsize = 7, col = "black"))
}

### Basic Rank plot ####
plot(f.results, y = txpRanks(f.results), 
     labels = c(1:10), 
     pch = 1, size = unit(0.75, 'char') #data point shape and size
)

### Data preparation for Upgrade Rank plot ####
library(dplyr)
# Load selected chemical info

Liver_Tox_1 <- Liver_Tox %>%
  mutate(id = row_number()) %>%
  select(id, DTXSID, PREFERRED_NAME.x)

# Transform matrix to dataframe
toxpi_scores <- as.data.frame(toxpi_scores) 

# Add toxpi_score variable which is the sum of all slices
toxpi_scores <- toxpi_scores %>%
  mutate(id = row_number()) %>%
  rowwise() %>%
  mutate(toxpi_score = sum(c_across(Slice1:Slice14))) %>%
  ungroup()

# Merged chemcial info and slice, toxpi scores
Liver_Toxpi_Slices <- Liver_Tox_1 %>%
  full_join(toxpi_scores, by = 'id') %>%
  rename(PREFERRED_NAME = PREFERRED_NAME.x)

# Output ToxpiScore
library(writexl)
write_xlsx(Liver_Toxpi_Slices,'8-Liver_Tox_Slices_N=94.xlsx')

#### Upgrade Rank plot ####
library(ggplot2)

RankPlot <- ggplot(Liver_Toxpi_Slices, 
                   aes(x = toxpi_score, y = reorder(PREFERRED_NAME, toxpi_score))) +
  geom_point() +
  scale_x_continuous(limits = c(0.2, 0.75), breaks = seq(0.2, 0.75, 0.05)) +
  labs(x = "Tox Score", y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(
      hjust = 1,  
      size = 6.0    
    )
  )

show(RankPlot)

# Output plot
ggsave("8.2_RankPlot_Tox_Liver.pdf", plot = RankPlot, 
       width = 6, height = 8, units = "in", dpi = 600)

### Hierarchical Clustering ####
f.hc <- hclust(dist(txpSliceScores(f.results)))
plot(f.hc, hang = -1, labels = txpIDs(f.results), xlab="", sub="")


####### 9.ExpTox_Liver_BEQ_BER_RI_Calculation----------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

#### API-CompTox-ctxR package, extract Liver_Xenobiotic metabolism data ####

library(ctxR)
### API key will be stored as the variable my_key
my_key <- '5e944183-665f-4e1b-9ff0-c7e27b8b09e4'

### Retrieve all assays
assays <- get_all_assays(API_key = my_key) #1499 assays

## Select human, liver, nuclear receptor
library(dplyr)
assays_1 <- assays %>%
  filter(organism == 'human' &
           tissue == 'liver' &
           intendedTargetFamilySub == 'xenobiotic metabolism') #34 assays

XM_aeid <- assays_1 %>%
  select(aeid)

library(writexl)
write_xlsx(assays_1, '9.1-Liver_XenoMetabolism_assays.xlsx')

# Extract DTXSID from Liver 
Liver <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                    sheet = 1)
DTXSID <- Liver %>%
  select(DTXSID) %>%
  distinct() %>%
  pull(DTXSID) #94 DTXSID

### Retrieve bioactivity data from NR_AEID batch
XM_bio_Liver <- get_bioactivity_details_batch(
  API_key = my_key,
  DTXSID = DTXSID #Liver compounds
)


### Merge multiple lists into a database
library(dplyr)
library(purrr)

# Function to check if a dataframe is valid (non-empty)
is_valid_df <- function(df) {
  if (is.data.frame(df) && nrow(df) > 0 && ncol(df) > 0) {
    return(TRUE)
  }
  return(FALSE)
}

# Filter out empty dataframes from the list
XM_bio_Liver_filtered <- map(XM_bio_Liver, ~ {
  if (is_valid_df(.x)) .x else NULL
})

# Remove NULL elements from the list
XM_bio_Liver_filtered <- compact(XM_bio_Liver_filtered)

XM_bio_Liver_merge <- XM_bio_Liver_filtered %>%
  map_df(~ .x) #49215 records

# XM_aeid_bio_Liver
XM_bio_Liver_merge_1 <- XM_bio_Liver_merge %>%
  select(dtxsid, aeid, ac50)

XM_aeid_bio_Liver <- left_join(XM_aeid, XM_bio_Liver_merge_1, by = 'aeid')


unique(XM_aeid_bio_Liver$dtxsid) #35 compounds

unique(XM_aeid_bio_Liver$aeid) #34 XM assays

# Save data
library(writexl)
write_xlsx(XM_aeid_bio_Liver, '9.2-XM_aeid_bio_Liver_N=35.xlsx') 

#### BEQ calcualtion ####
library(readxl)
Liver_XM_bio <- read_excel('9.2-XM_aeid_bio_Liver_N=35.xlsx',
                           sheet = 1) #548 records


## calculate ac50_ref for assay: select min ac50 -----------

ac50_ref_df <- Liver_XM_bio %>%
  group_by(aeid) %>%
  dplyr::summarise(
    # 
    ac50_ref_min = min(ac50, na.rm = TRUE)) %>%
  ungroup()

## load CB data ----------------------------------------------------
load('5-Liver_Pred_Expo_CytoAC50_N=94.RData')

Liver_CB <- Liver_Pred_Expo_CytoAC50 %>%
  select(DTXSID, CB_P50_uM) %>%
  distinct() %>%
  filter(!is.na(CB_P50_uM)) %>%
  dplyr::rename(dtxsid = DTXSID)

# merge CB and ac50_ref
df_base <- Liver_XM_bio %>%
  left_join(Liver_CB,   by = "dtxsid") %>%
  left_join(ac50_ref_df, by = "aeid")


## Calculate CB/AC50, BEQ, BEQ% -----------------------------------
df_beq <- df_base %>%
  mutate(CB_AC50 = CB_P50_uM / ac50) %>%
  mutate(BEQ_min = CB_AC50 * ac50_ref_min) %>%
  group_by(aeid) %>%
  mutate(sum_BEQ_min = sum(BEQ_min, na.rm = TRUE),
         BEQ_percent_min = BEQ_min / sum_BEQ_min * 100) %>%
  ungroup()

# Each compound conresponding to multiple target assays
library(tidyr)
df_beq_1 <- df_beq %>%
  select(dtxsid, aeid, BEQ_percent_min) %>%
  pivot_wider(names_from = aeid,
              values_from = BEQ_percent_min) 

# select assays 
df_beq_2 <- df_beq_1 %>%
  select(1:2, any_of(as.character(962:988))) #15 assays

Liver_XM_BEQ_CB <- df_beq_2
save(Liver_XM_BEQ_CB, file = '9.1-Liver_XM_BEQ_CB_N=35.RData')


### RI calculation using SEEM3 ####
library(readxl)
Toxpi <- read_excel('8-Liver_Tox_Slices_N=94.xlsx',
                    sheet = 1)

load('5-Liver_Pred_Expo_CytoAC50_N=94.RData')

library(dplyr)
Toxpi <- Toxpi %>%
  select(DTXSID, PREFERRED_NAME, toxpi_score)

SEEM3 <- Liver_Pred_Expo_CytoAC50 %>%
  select(DTXSID, Ex_P50) %>%
  distinct() %>%
  filter(!is.na(Ex_P50)) %>%
  mutate(log_seem3 = log(Ex_P50)) # log-transformed SEEM3

df <- left_join(Toxpi, SEEM3, by = 'DTXSID')

## Normalized SEEM3

df$seem3_normalized <- (df$log_seem3 - min(df$log_seem3, na.rm = TRUE)) / 
  (max(df$log_seem3, na.rm = TRUE) - min(df$log_seem3, na.rm = TRUE))

## Calulate RI
df <- df %>%
  mutate(RI = seem3_normalized * toxpi_score) #87 compounds

# Save data
Liver_RI <- df
save(Liver_RI, file = '9.2-Liver_RI_N=87.RData')

##### BER Calculation using HTTK #####
library(readxl)

Che_list <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                       sheet = 1)

### HTTK model to calculate Css ###
library(httk)
library(readxl)
library(ggplot2)
library(dplyr)

### Step 1: Loading In Vitro Screening Data ####

load('invitrodb_3_5_mc5.Rdata')

# Built new subset of chemicals with human HTTK data
toxcast.httk <- subset(mc5, dsstox_substance_id %in% get_cheminfo(
  info='DTXSID', suppress.messages = TRUE))

# Extract 'DTXSID' from our chemcial list
my.chems <- Che_list$DTXSID #76

# Extract cooresponding chemicals from Toxcast
example.toxcast <- as.data.frame(mc5[mc5$dsstox_substance_id %in% my.chems,]) #51074 records

# pick 9 variables
example.toxcast <- example.toxcast[, c("chnm",
                                       "dsstox_substance_id",
                                       "spid",
                                       "hitc",
                                       "modl",
                                       "aeid",
                                       "modl_ga",
                                       "modl_ac10",
                                       "modl_acc")]

# reduce precision to decrease space:
cols <- c("modl_ga", "modl_ac10", "modl_acc")
for (this.col in cols)
  example.toxcast[, this.col] = signif(example.toxcast[, this.col], 3)

#Each chemical might have multiple assays for which there is a hit. 
#Each chemical may have been run multiple times (experimental replicates), 
#so that there may be multiple samples from the same chemical for the same assay. 
#Different samples are indicated by the column “spid” (sample id).

knitr::kable(head(example.toxcast), caption = "ToxCast In Vitro Bioactivity Data",
             floating.environment="sidewaystable")

# Find the in vitro concentrations for compounds
toxcast.table <- NULL 
for (this.id in unique(example.toxcast$dsstox_substance_id))
{
  this.subset <- subset(example.toxcast, dsstox_substance_id == this.id)
  these.hits <- subset(this.subset, hitc==1)
  if (dim(these.hits)[1]>0){
    this.row <- data.frame(Compound=as.character(this.subset[1,"chnm"]),
                           DTXSID=this.id,
                           Total.Assays = dim(this.subset)[1],
                           Unique.Assays = length(unique(this.subset$aeid)),
                           Total.Hits = dim(these.hits)[1],
                           Unique.Hits = length(unique(these.hits$aeid)),
                           Low.AC50 = signif(min(these.hits$modl_ga),3),
                           Low.AC10 = signif(min(these.hits$modl_ac10),3),
                           Low.ACC = signif(min(these.hits$modl_acc),3),
                           Q10.AC50 = signif(quantile(these.hits$modl_ga,probs=0.1),3),
                           Q10.AC10 = signif(quantile(these.hits$modl_ac10,probs=0.1),3),
                           Q10.ACC = signif(quantile(these.hits$modl_acc,probs=0.1),3),
                           Med.AC50 = signif(median(these.hits$modl_ga),3),
                           Med.AC10 = signif(median(these.hits$modl_ac10),3),
                           Med.ACC = signif(median(these.hits$modl_acc),3),
                           Q90.AC50 = signif(quantile(these.hits$modl_ga,probs=0.9),3),
                           Q90.AC10 = signif(quantile(these.hits$modl_ac10,probs=0.9),3),
                           Q90.ACC = signif(quantile(these.hits$modl_acc,probs=0.9),3)
    )
    toxcast.table <- rbind(toxcast.table, this.row)
  }
}
rownames(toxcast.table) <- seq(1,dim(toxcast.table)[1]) #39
knitr::kable(head(toxcast.table[,1:6]), caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")

### Step 2: HTTK deduce Css_Plasma with in vitro ####

#Steady State Plasma Concentration Css 
#calc_mc_css(): predict Css resulting from an ongoing 1 mg/kg/day exposure.

# Assume toxcast.table already has columns "Css_95th", "Css_50th", "Css_5th"
for (this.id in unique(toxcast.table$DTXSID)) {
  # Check if a chemical substance is in the HTTK database
  if (this.id %in% get_cheminfo(info = "dtxsid", suppress.messages = TRUE)) {
    # 
    set.seed(12345)
    
    # Calculating Css at different percentiles
    css_results <- calc_mc_css(dtxsid = this.id,
                               output.units = "uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # 
                               suppress.messages = TRUE)
    
    # Store the results in a table
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_5th"]  <- css_results[3]
    
    # Tag Css Type
    toxcast.table[toxcast.table$DTXSID == this.id, "Css.Type"] <- "in vitro"
  }
}

#Let’s look at just the relevant columns from our table. 
#Any Css values of “NA” (not a number) indicate that we don’t currently have in vitro TK data for those chemicals:

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th', "Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")


### Step 3: HTTK deduce Css_Plasma with QSPR Predictions ####

#For chemicals where no in vitro toxicokinetic data are available, 
#we can sometimes load predictions from in silico models. 

load_sipes2017()

library(httk)

# 
if (!"Css_95th" %in% colnames(toxcast.table)) {
  toxcast.table$Css_95th <- NA
}
if (!"Css_50th" %in% colnames(toxcast.table)) {
  toxcast.table$Css_50th <- NA
}
if (!"Css_5th" %in% colnames(toxcast.table)) {
  toxcast.table$Css_5th <- NA
}
if (!"Css.Type" %in% colnames(toxcast.table)) {
  toxcast.table$Css.Type <- NA
}

# Iterate over each chemical ID
for (this.id in unique(toxcast.table$DTXSID)) {
  
  # Check if the chemical is in the HTTK database and Css is empty
  if (this.id %in% get_cheminfo(info="dtxsid", suppress.messages=TRUE) &
      is.na(toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"])) {
    
    # Set the random number seed to ensure repeatable results
    set.seed(12345)
    
    # Calculate 95th, 50th and 5th percentile Css
    css_results <- calc_mc_css(dtxsid=this.id,
                               output.units="uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # Default number of samples
                               suppress.messages=TRUE)
    
    # Store calculation results in different columns
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_5th"] <- css_results[3]
    
    # Mark Css type as "in silico"
    toxcast.table[toxcast.table$DTXSID==this.id, "Css.Type"] <- "in silico"
  }
}

#Take another look at our table – many of the gaps are now filled in.

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th',"Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")

sum(!is.na(toxcast.table$Css_50th))

## Output
Liver_Css_Plasma = toxcast.table %>%
  filter(!is.na(Css_50th)) #47 compounds

save(Liver_Css_Plasma, file = '9.3-Liver_Css_Plasma_N=47.RData')

## Step 4: Reverse Dosimetry In Vitro-In Vivo Extropolation ####
load('9.3-Liver_Css_Plasama_N=47.RData')
#Don’t forget the ToxCast concentrations are given on the log-10 scale:
Liver_Css_Plasma$EquivDose <- signif(10^Liver_Css_Plasma$Q10.AC50 / Liver_Css_Plasma$Css_50th,
                                     3)

knitr::kable(Liver_Css_Plasma[1:10,c("Compound","Q10.AC50","Css_50th","EquivDose")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")          

## Step 5: Compare with Estimate Exposure Rates ####

# Load SEEM3 database
load('/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/External_Prediction.RData')

library(dplyr)
SEEM <- External_Prediction %>%
  select(DTXSID, InChIKey, P5, P50, P95)

# Pick interested chemicals
example.seem <- as.data.frame(subset(SEEM,
                                     DTXSID %in% Liver_Css_Plasma$DTXSID))

# Merge database
Liver_Css_Plasma <- left_join(Liver_Css_Plasma, example.seem, by = 'DTXSID')

#Calculate a bioactivity:exposure ratio (BER)
Liver_Css_Plasma$BER <- signif(Liver_Css_Plasma$EquivDose/Liver_Css_Plasma$P95,3)

Liver_BER <- Liver_Css_Plasma %>%
  select(InChIKey, DTXSID, BER) %>%
  filter(!is.na(BER))

unique(Liver_BER$InChIKey)
unique(Liver_BER$DTXSID) #47 compounds

save(Liver_BER, file = '9.4-Liver_HTTK_BER_N=47.RData')




####### 10.ExpTox_Plot_Liver----------------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

### Load data
library(dplyr)

# RI data
load('9.2-Liver_RI_N=87.RData')

# BEQ data
load('9.1-Liver_XM_BEQ_CB_N=35.RData')

BEQ <- Liver_XM_BEQ_CB %>%
  rename(DTXSID = dtxsid) %>%
  rename_with(~ paste0("BEQ_", .), .cols = -1) #35 compounds

# BER data
load('9.4-Liver_HTTK_BER_N=47.RData')
BER <- Liver_BER %>%
  select(DTXSID, BER) #47 compounds

#### Merge RI 和 BEQ, BER ####
df <- Liver_RI %>%
  left_join(BER, by = 'DTXSID') %>%
  left_join(BEQ, by = 'DTXSID') 

unique(df$RI)
unique(df$BER)

#### ExpTox calculation ####
library(toxpiR)

# View all functions
lsf.str("package:toxpiR")

## Create slice with transformation functions
# Slice 1: BEQ, # log10(x)-log10(min)
# Slice 2: BER, # -log10(x)+log10(max)
# Slice 3: RI, # linear

# Automatically extract all column names starting with 'BEQ_'
BEQ_columns <- grep("^BEQ_", names(df), value = TRUE)

f.slices <- TxpSliceList(
  Slice1 = TxpSlice(BEQ_columns),
  Slice2 = TxpSlice('BER'),
  Slice3 = TxpSlice('RI')
)

## Method#1: Not consider missing value
linear <- function(x) {
  # Example linear transformation
  return(x)
}

final.trans <- TxpTransFuncList(
  f1 = function(x) {
    if (all(is.na(x))) return(x)
    log10(x) - log10(min(x, na.rm = TRUE))
  },
  f2 = function(x) {
    if (all(is.na(x))) return(x)
    -log10(x) + log10(max(x, na.rm = TRUE))
  },
  f3 = linear
)


## Goal - Create ToxPi model.
# All Slices, weight = 1.
f.model <- TxpModel(txpSlices = f.slices,
                    txpWeights = c (1,1,1),
                    txpTransFuncs = final.trans)


## Calculate ToxPi scores
f.results <- txpCalculateScores(model = f.model,
                                input = df,
                                id.var = 'PREFERRED_NAME')
# ToxPi scores
toxpi_scores <- txpSliceScores(f.results)
print(toxpi_scores)

#### Pie Plotting ####

colors <- c("#E69F00", "#56B4E9", "#009E73")

library(grid)
plot(f.results,
     showScore = TRUE,
     fills = colors # fill colors
)

#### Edit Pie Plotting ####
# List all grobs in the current graph
grid.ls()
# Find all grobs containing 'label' and change the font size
all_grobs <- grid.ls(print = FALSE)


## Modify the font size of all labels
label_grobs <- grep("label", all_grobs$name, value = TRUE)

for (label in label_grobs) {
  grid.edit(label, gp = gpar(fontsize = 6))
}

## Edit 'Toxpi Scores' characteristics
for (s in sprintf("pie-%d-radSum", 1:94)) {
  grid.edit(s, gp = gpar(fontsize = 7, col = "black"))
}


### Basic Rank plot ####
plot(f.results, y = txpRanks(f.results), 
     labels = c(1:10), 
     pch = 1, size = unit(0.75, 'char') #data point shape and size
)

### Data preparation for Upgrade Rank plot ####
library(dplyr)
# Load selected chemical info
df_1 <- df %>%
  mutate(id = row_number())

df_2 <- df_1 %>%
  mutate(id = row_number()) %>%
  select(id, DTXSID, PREFERRED_NAME)

# Transform matrix to dataframe
toxpi_scores <- as.data.frame(toxpi_scores) 

# Add toxpi_score variable which is the sum of all slices
toxpi_scores <- toxpi_scores %>%
  mutate(id = row_number()) %>%
  rowwise() %>%
  mutate(toxpi_score = sum(c_across(Slice1:Slice3))) %>%
  ungroup()

# Merged chemical info and slice, toxpi scores
Liver_ExpTox_Slices <- df_2 %>%
  full_join(toxpi_scores, by = 'id')

library(writexl)
write_xlsx(Liver_ExpTox_Slices, '10-Liver_ExpTox_Slice_N=94.xlsx')

#### Upgrade Rank plot ####
library(ggplot2)

# 
Liver_ExpTox_Slices <- Liver_ExpTox_Slices %>%
  filter(!is.na(toxpi_score))

RankPlot <- ggplot(Liver_ExpTox_Slices, 
                   aes(x = toxpi_score, y = reorder(PREFERRED_NAME, toxpi_score))) +
  geom_point() +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.1)) +
  labs(x = "ExpTox Score", y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(
      hjust = 1,  # 
      size = 6.0    # 
    )
  )

show(RankPlot)

# Output plot
ggsave("10.2-RankPlot_ExpTox_Liver.pdf", plot = RankPlot, 
       width = 6, height = 8, units = "in", dpi = 600)

### Hierarchical Clustering ####
f.hc <- hclust(dist(txpSliceScores(f.results)))
plot(f.hc, hang = -1, labels = txpIDs(f.results), xlab="", sub="")

####### 11.Plot_Rank_Change_Liver-----------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')
library(readxl)
library(dplyr)

Tox <- read_excel('8-Liver_Tox_Slices_N=94.xlsx',
                  sheet = 1)
Tox_1 <- Tox %>%
  select(PREFERRED_NAME, toxpi_score) %>%
  mutate(rank_tox = rank(-toxpi_score, ties.method = "average")) %>%
  rename(Tox_score = toxpi_score)

ExpTox <- read_excel('10-Liver_ExpTox_Slice_N=94.xlsx',
                     sheet = 1)
ExpTox_1 <- ExpTox %>%
  select(PREFERRED_NAME, toxpi_score) %>%
  mutate(rank_exp = rank(-toxpi_score, ties.method = 'average')) %>%
  rename(ExpTox_score = toxpi_score)

Tox_Exp <- left_join(Tox_1, ExpTox_1, by = 'PREFERRED_NAME')

Tox_Exp <- Tox_Exp %>%
  mutate(change = rank_tox - rank_exp)

# Plotting ranking changes

library(ggplot2)
library(ggrepel)

ggplot(Tox_Exp, aes(x = rank_tox, y = rank_exp, label = PREFERRED_NAME)) +
  geom_point(aes(color = change), size = 1.5) +
  geom_text_repel(size = 3.5, max.overlaps = Inf) + #  geom_text_repel
  scale_color_gradient2(low = "blue", high = "red", mid = "purple", midpoint = 0) +
  labs(title = "Rank Change Plot", x = "Tox_rank", y = "ExpTox_rank") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),       
    axis.title = element_text(size = 12),                    
    axis.text = element_text(size = 12),                    
    legend.title = element_text(size = 11),                   
    legend.text = element_text(size = 11)                     
  )



