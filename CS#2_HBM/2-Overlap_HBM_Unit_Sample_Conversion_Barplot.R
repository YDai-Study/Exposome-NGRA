
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

## Load All_HBM_Merge data

load('/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/All_HBM_Merge_New.RData')

# Filter EU_HBM, Remove Non-EU Region (Israel)
library(dplyr)
EU_HBM <- All_HBM_Merge_New %>%
  filter(Region %in% c("Northern EU", "Western EU", "Southern EU", "Eastern EU"))

# Filter US_HBM
US_HBM <- All_HBM_Merge_New %>%
  filter(Country == "United States")

# Filter Canada_HBM
CA_HBM <- All_HBM_Merge_New %>%
  filter(Country == "Canada")

## Load overlapping 76 compounds
library(readxl)
Overlap_HBM <- read_excel('1-Overlap_EU_US_CA_N=76.xlsx', 
                          sheet = 1)
Overlap_HBM_1 <- Overlap_HBM %>%
  select(InChIKey, Chemical_Group)

##### 1. Extract 76 compounds concentration from HBM data  ####

# Extract INCHIKEY List
INC_List <- Overlap_HBM$InChIKey[1:76]

# EU region
EU_HBM_1 <- EU_HBM[EU_HBM$InChIKey %in% INC_List,]

# US region
US_HBM_1 <- US_HBM[US_HBM$InChIKey %in% INC_List,]

# Canada
CA_HBM_1 <- CA_HBM[CA_HBM$InChIKey %in% INC_List,]

## Merge
EU_US_CA_HBM <- bind_rows(EU_HBM_1, US_HBM_1, CA_HBM_1) %>%
  select(-Chemical_Group)

## Input Chemical_Group
Overlap_HBM_Conc <- left_join(EU_US_CA_HBM, Overlap_HBM_1,
                                by = 'InChIKey')

library(writexl)
write_xlsx(Overlap_HBM_Conc, '2.1-Overlap_HBM_Conc.xlsx')

##### 2. Select concentration statistic values ####
library(dplyr)
Overlap_HBM_Conc_1 <- Overlap_HBM_Conc %>%
  mutate(Concentration = coalesce(P50, `Geometric Mean`, Mean, P75, P90, P95)) %>%
  select(-P05, -P10, -P25, -P50, -P75, -P90, -P95, 
         -Minimum, Mean, -`Geometric Mean`, -Maximum) %>%
  relocate(Chemical_Name, InChIKey, CAS_Number, DTXSID, Average_Mass,
           Chemical_Name_Original, Chemical_Group, Concentration, everything()) %>%
  dplyr::rename(Unit = Concentration_Unit)

##### 3. Concentration unit conversion ####

unique(Overlap_HBM_Conc_1$Unit)
colnames(Overlap_HBM_Conc_1)

# Unit: "ng/mL" 
# Unit: "ng/L", "pg/mL" change to  "µg/L", Conc/1000
# Unit: "ng/g" change to  "µg/g", Conc/1000
# Unit: "ng/g creatinine" change to "µg/g crt", Conc/1000

# Unit: "µg/dL" change to "µg/L", Conc * 10
# Unit: "µg/g lipid" change to "ng/g lipid", Conc*1000

Overlap_HBM_Conc_2 <- Overlap_HBM_Conc_1 %>%
  mutate(Unit = ifelse(Unit == "ng/mL", "µg/L", Unit),
         
         Concentration = ifelse(Unit == "ng/L", Concentration / 1000, Concentration),
         Unit = ifelse(Unit == "ng/L", "µg/L", Unit),
         
         Concentration = ifelse(Unit == "pg/mL", Concentration / 1000, Concentration),
         Unit = ifelse(Unit == "pg/mL", "µg/L", Unit),
         
         Concentration = ifelse(Unit == "ng/g", Concentration / 1000, Concentration),
         Unit = ifelse(Unit == "ng/g", "µg/g", Unit),
         
         Concentration = ifelse(Unit == "ng/g creatinine", Concentration / 1000, Concentration),
         Unit = ifelse(Unit == "ng/g creatinine", "µg/g creatinine", Unit),
         
         Concentration = ifelse(Unit == "µg/dL", Concentration * 10, Concentration),
         Unit = ifelse(Unit == "µg/dL" , "µg/L", Unit),
         
         Concentration = ifelse(Unit == "µg/g lipid", Concentration * 1000, Concentration),
         Unit = ifelse(Unit == "µg/g lipid", "ng/g lipid", Unit)
         ) %>%
  mutate(Unit = ifelse(Unit == "µg/g sample", "µg/g breast milk", Unit))

unique(Overlap_HBM_Conc_2$Unit)

##### 4. Sample type ####
unique(Overlap_HBM_Conc_2$Sample_Type)
unique(Overlap_HBM_Conc_2$Sample)

Overlap_HBM_Conc_3 <- Overlap_HBM_Conc_2 %>%
  mutate(
    Sample_INF = case_when(
      Sample == "Human (Urine)" ~ "Urine",
      Sample == "urine" ~ "Urine",
      
      Sample == "Human (Hair)" ~ "Hair",
      
      Sample == "Human (Blood)" ~ "Blood",
      Sample == "Serum" ~ "Blood",
      Sample == "serum" ~ "Blood",
      Sample == "blood" ~ "Blood",
      Sample == "plasma" ~ "Blood",
      Sample == "pooled serum" ~ "Blood",
      
      Sample == "Human (Cord Blood)" ~ "Cord blood",
      Sample == "Placenta tissue (ng/g)" ~ "Placenta",
      Sample == "Breast Milk" ~ "Breast milk",
      Sample_Type %in% c(
        "cord whole blood", "cord plasma", "cord plasma fat",
        "Cord Blood Plasma", "Cord Blood Whole Blood", "Cord Blood Serum"
      ) ~ "Cord blood",
      TRUE ~ Sample  # 默认保留 Sample 原始值
    )
  )

unique(Overlap_HBM_Conc_3$Sample_INF)
unique(Overlap_HBM_Conc_3$Sample_Type)

## Save data
library(writexl)
write_xlsx(Overlap_HBM_Conc_3, '2.2-Overlap_HBM_Conc_Unit_Conversion.xlsx')

##### 5. Barplot_Overlap_Sample ####
library(readxl)
HBM <- read_excel('2.2-Overlap_HBM_Conc_Unit_Conversion.xlsx',
                  sheet = 1)
# Delete repeated records for same values in INCHIKEY and Sample
HBM_unique <- HBM %>% 
  distinct(InChIKey, Sample_INF) #146

# Transform chr to factor
HBM_unique$Sample_INF <- as.factor(HBM_unique$Sample_INF)
unique(HBM_unique$Sample_INF)

## Barplot: compound counts in different samples 

# calculate the counts of inchikey for each sample
library(dplyr)

INCHIKEY_counts <- HBM_unique %>%
  group_by(Sample_INF) %>%
  summarise(INCHIKEY_counts = n_distinct(InChIKey))

library(ggplot2)

ggplot(INCHIKEY_counts, aes(x = reorder(Sample_INF, -INCHIKEY_counts), 
                            y = INCHIKEY_counts, fill = Sample_INF)) +
  geom_bar(stat = "identity", width = 0.8) +  
  geom_text(aes(label = INCHIKEY_counts), hjust = -0.3, size = 5) +  # Label the bars with values
  theme_minimal() +
  labs(x = "", y = "Count of compounds by sample") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12),
        legend.position = "none") +  
  scale_y_discrete(position = "left") +  
  scale_y_reverse() + 
  coord_flip()  # Flip the axes so the bars appear horizontally


