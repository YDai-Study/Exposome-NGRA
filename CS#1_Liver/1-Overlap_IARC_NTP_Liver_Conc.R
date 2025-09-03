
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

##### 1. Overlap Liver compounds between NTP and IARC list ####
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

##### 2. Extract 94 compounds from Exposome data #####

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

##### 3. Check 94 compounds in Exposome data #####

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
