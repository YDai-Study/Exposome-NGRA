
###### 1.Overlap_EU_US_CA+DonutChart--------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

# Load All_HBM_Merge data
load('/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/All_HBM_Merge_New.RData')
unique(All_HBM_Merge_New$Region)
unique(All_HBM_Merge_New$Country)

# Filter EU_HBM, Remove Non-EU Region (Israel)
library(dplyr)
EU_HBM <- All_HBM_Merge_New %>%
  filter(Region %in% c("Northern EU", "Western EU", "Southern EU", "Eastern EU"))

# Filter US_HBM
US_HBM <- All_HBM_Merge_New %>%
  filter(Country == "United States")

# Filter Canada_HBM
Canada_HBM <- All_HBM_Merge_New %>%
  filter(Country == "Canada")

##### 1.1.Overlap Results in EU,US,CA HBM Databases ######
# create list
EU_list <- EU_HBM$InChIKey[1:193472]
US_list <- US_HBM$InChIKey[1:28460]
CA_list <- Canada_HBM$InChIKey[1:12342]

#create a list include all list value
venn_list <- list (EU_HBM = EU_list,
                   CA_HBM = CA_list,
                   US_HBM = US_list) 

#remove NA value from every vector in list
venn_list = purrr::map(venn_list, na.omit)


## Overlapping between HBM Databases

library(ggVennDiagram)
library(ggplot2)

ggVennDiagram(venn_list[1:3], label = c('count'), 
              label_alpha = 0,
              label_size = 5,
              edge_size = 0.1,
              set_size = 5) + scale_fill_gradient(low="grey90",high = "#488AC7")

process_region_data(Venn(venn_list[1:3]))

ggVennDiagram(venn_list[1:3], force_upset = TRUE, order.set.by = "name", 
              order.intersect.by = "size")


## View intersection results
library(VennDiagram)
library(grid)
library(futile.logger)

inter <- get.venn.partitions(venn_list[1:3])
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
inter <- subset(inter, select = -..values.. )
inter <- subset(inter, select = -..set.. )

### Export intersection results ####
write.table (inter, '1-Overlap_EU_US_CA.csv', 
             row.names = FALSE,
             sep = ',',
             quote = FALSE)

#76 overlapping compounds
INC_76_string <- inter$values[1]
INC_76_vector <- unlist(strsplit(INC_76_string, "\\|"))
INC_76_df <- data.frame(InChIKey = trimws(INC_76_vector))

## Load CompTox Data, link compound identifiers
load('DSSTox_Feb_2024.RData')
library(dplyr)

INC_76_df <- DSSTox_Feb_2024[DSSTox_Feb_2024$INCHIKEY
                             %in% INC_76_vector, ]
# remove duplicated INCHIKEY
INC_76_df_1 <- INC_76_df[-c(72, 69), ] %>%
  dplyr::rename(InChIKey = INCHIKEY)


## To find Chemical_Group from Canada_HBM
INC_76_df_2 <- Canada_HBM[Canada_HBM$InChIKey %in% INC_76_vector,]


INC_76_df_2 <- INC_76_df_2 %>%
  select(InChIKey, Chemical_Group) %>%
  distinct()

INC_76_df_3 <- full_join(INC_76_df_1, INC_76_df_2, by = 'InChIKey')

# output 76 overlapping compounds
library(writexl)
write_xlsx(INC_76_df_3, '1-Overlap_EU_US_CA_N=76.xlsx')



##### 1.2. DonutChart_HBM_ChemGroup_N=76 ####
library(readxl)
HBM <- read_excel('/Users/ydai/Desktop/Exposome_NGRA/2505-HBM_Analysis/1-Overlap_EU_US_CA_N=76.xlsx',
                  sheet = 1)

## Variable transformation
HBM$Chemical_Group <- as.factor(HBM$Chemical_Group)

## Pie chart: percentage 
# Data transformation
library(dplyr)
library(scales)
library(ggplot2)

df <- HBM %>%
  group_by(Chemical_Group) %>%
  count() %>%
  ungroup() %>%
  mutate(perc = n / sum(n)) %>%
  arrange(desc(n)) %>%  # Sort by n variables in descending order
  mutate(labels = scales::percent(perc))

print(df)

## Shows percentage
ggplot(df, aes(x = "", y = perc, fill = Chemical_Group)) +
  geom_col() +
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") 

## Shows Number
# Generate color vector
cols <- hcl.colors(length(unique(HBM$Chemical_Group)), "Temps")

# Use ggplot2 to draw a hollow pie chart, sorting Chemical_Group in descending order by n
ggplot(df, aes(x = 2, y = n, fill = reorder(Chemical_Group, -n))) +  # Sort by n in descending order Chemical_Group
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = cols) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
  xlim(0.5, 2.5) +  # Adjust the x-axis range to create a hollow effect
  theme(legend.position = "right") +
  labs(fill = "Chemical Group", title = "Chemical Group Distribution")







###### 2.Overlap_HBM_Unit_Sample_Conversion-------------------------------------

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

##### 2.1. Extract 76 compounds concentration from HBM data  ####

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

##### 2.2. Select concentration statistic values ####
library(dplyr)
Overlap_HBM_Conc_1 <- Overlap_HBM_Conc %>%
  mutate(Concentration = coalesce(P50, `Geometric Mean`, Mean, P75, P90, P95)) %>%
  select(-P05, -P10, -P25, -P50, -P75, -P90, -P95, 
         -Minimum, Mean, -`Geometric Mean`, -Maximum) %>%
  relocate(Chemical_Name, InChIKey, CAS_Number, DTXSID, Average_Mass,
           Chemical_Name_Original, Chemical_Group, Concentration, everything()) %>%
  dplyr::rename(Unit = Concentration_Unit)

##### 2.3. Concentration unit conversion ####

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

##### 2.4. Sample type ####
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
      TRUE ~ Sample  
    )
  )

unique(Overlap_HBM_Conc_3$Sample_INF)
unique(Overlap_HBM_Conc_3$Sample_Type)

## Save data
library(writexl)
write_xlsx(Overlap_HBM_Conc_3, '2.2-Overlap_HBM_Conc_Unit_Conversion.xlsx')

##### 2.5. Barplot_Overlap_Sample ####
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





####### 3.Overlap_HBM_SankeyDiagram---------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

## Load HBM overlapping compounds with selected concentration
library(readxl)
HBM <- read_excel('2.2-Overlap_HBM_Conc_Unit_Conversion.xlsx',
                  sheet = 1)
unique(HBM$InChIKey) #76 INCHIKEYs

library(dplyr)
# Input a new variable
HBM <- HBM %>%
  mutate(Conc_LOD = case_when(
    is.na(Concentration) ~ '<LOD',
    !is.na(Concentration) ~ '>LOD'))

# Mark records with conc_lod value <lod
HBM <- HBM %>%
  mutate(conc_lod_flag = ifelse(Conc_LOD == "<LOD", 1, 0))

# Group by inchikey and sample, and prioritize deleting records with a conc_lod_flag value of 1.
HBM_cleaned <- HBM %>%
  group_by(InChIKey, Sample_INF) %>%
  arrange(conc_lod_flag) %>% # Sort by conc_lod_flag first, make sure <lod is in front
  slice(1) %>% # Keep the first record of each group
  ungroup() %>% # Ungroup
  select(-conc_lod_flag) 

unique(HBM_cleaned$InChIKey) #76 INCHIKEYs

##### Check the data ####
# show the INCHIKEY with Conc <LOD
data <- HBM_cleaned %>%
  filter(is.na(Concentration))
unique(data$InChIKey) #10 INCHIKEYs

# show the INCHIKEY with Conc >LOD
data.1 <- HBM_cleaned %>%
  filter(Concentration != 'NA')
unique(data.1$InChIKey) #76 INCHIKEYs

#### Overlap: <LOD and >LOD ####
# Overlap between data and data.1
list <- data$InChIKey[1:17]
list.1 <- data.1$InChIKey[1:129]

venn_list <- list(Below_LOD = list,
                  Above_LOD = list.1)
#remove NA value from every vector in list
venn_list = purrr::map(venn_list, na.omit)

library(ggvenn)
library(VennDiagram)
library(futile.logger)

## 
ggvenn(
  venn_list, 
  fill_color = c("#EFC000FF", "#0073C2FF"),text_size = 5,
  stroke_size = 0.1, set_name_size = 5, show_percentage = FALSE)

inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '|')
inter <- subset(inter, select = -..values.. )
inter <- subset(inter, select = -..set.. )

#export intersection results
write.table (inter, '3-Overlap_HBM_Below_Above_LOD.csv', 
             row.names = FALSE,
             sep = ',',
             quote = FALSE)


##### Sankey Chart_Original ####
#install.packages("remotes")
#remotes::install_github("davidsjoberg/ggsankey")

library(ggsankey)
library(ggplot2)
library(dplyr)

# Step 1
df <- HBM_cleaned %>%
  make_long(Sample_INF,Chemical_Group, Conc_LOD)
df

### Sankey Chart: show the count of each node 

# Step 2
dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()


# Step 3
df2 <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

# Chart 2
pl <- ggplot(df2, aes(x = x
                      , next_x = next_x
                      , node = node
                      , next_node = next_node
                      , fill = factor(node)
                      
                      , label = paste0(node," n=", n)
)
) 
pl <- pl +geom_sankey(flow.alpha = 0.5, show.legend = TRUE)
pl <- pl +geom_sankey_label(size = 3, color = "white", fill= "gray40", 
                            hjust = -0.1, nudge_y = 0.1)

pl <- pl +  theme_bw()
pl <- pl + theme(legend.position = "none")
pl <- pl +  theme(axis.title = element_blank()
                  , axis.text.y = element_blank()
                  , axis.ticks = element_blank()  
                  , panel.grid = element_blank())
pl <- pl + scale_fill_viridis_d(option = "inferno")
pl <- pl + labs(title = "Sankey diagram using ggplot")
pl <- pl + labs(subtitle = "using  David Sjoberg's ggsankey package")
pl <- pl + labs(caption = "@techanswers88")
pl <- pl + labs(fill = 'Nodes')


pl


#### Sankey diagram_update ####
library(dplyr)
library(ggplot2)
library(ggsankey)
library(stringr)

# 1) Label the three columns: Left (Sample), Middle (Chemical Group), Right (LOD)
df2 <- df2 %>%
  mutate(group = case_when(
    x == "Sample_INF"     ~ node,   
    x == "Chemical_Group" ~ node,   
    x == "Conc_LOD"       ~ node,   
    TRUE                  ~ "other"
  ))

# 2) Middle column color (from the Temps palette)
chem_groups <- df2 %>% filter(x == "Chemical_Group") %>% pull(node) %>% unique()
cols_mid <- hcl.colors(length(chem_groups), "Temps")
names(cols_mid) <- chem_groups

# 3) Specified colors for left and right columns 
cols_left <- c(
  "Urine"       = "#EE82EE",
  "Placenta"    = "#1f77b4",
  "Hair"        = "#2ca",
  "Cord blood"  = "#3EA055",
  "Breast milk" = "#8B8000",
  "Blood"       = "#F67280"
)
cols_right <- c(
  ">LOD" = "#E1D9D1",
  "<LOD"  = "#FFF9E3"
)

# 4) Merge colors and fill in a default color for any missing groups
all_colors <- c(cols_left, cols_right, cols_mid)
missing <- setdiff(unique(df2$group), names(all_colors))
all_colors[missing] <- "grey70"

# 5) Drawing
pl <- ggplot(
  df2,
  aes(x = x, next_x = next_x, node = node, next_node = next_node,
      fill = group, label = paste0(node, " n=", n))
) +
  geom_sankey(flow.alpha = 0.5, show.legend = FALSE) +
  geom_sankey_label(size = 3, color = "white", fill = "gray50",
                    hjust = -0.1, nudge_y = 0.1) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = all_colors) +
  labs(
    fill = "Nodes"
  )

print(pl)



####### 4.HBM_Cyto_AC50---------------------------------------------------------

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

library(tidyr)
library(dplyr)

HBM_Cyto_AC50_unique <- HBM_Cyto_AC50 %>%
  distinct(dtxsid, aeid, ac50, .keep_all = TRUE) #3676 records

HBM_Cyto_AC50_unique <- HBM_Cyto_AC50_unique %>%
  group_by(dtxsid, aeid, spid) %>%
  mutate(record_id = row_number()) %>%
  ungroup()

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





####### 5.HBM_HTTK_HExp_SEEM3---------------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

## Load HBM Overlapping compounds N=76 ##
library(readxl)
Che_list <- read_excel('1-Overlap_EU_US_CA_N=76.xlsx',
                       sheet = 1)

### HTTK model to calculate Css ###
library(httk)
library(readxl)
library(ggplot2)

### Step 1: Loading In Vitro Screening Data ####

load('/Users/ydai/Library/CloudStorage/OneDrive-Personal/2-2023-UMCU-Postdoc/1-UMCU-Postdoc Project/4-Prediction Tools/HTTK_R_Code/invitrodb_3_5_mc5.Rdata')

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
write_xlsx(toxcast.table, '5.1-HBM_HTTK_toxcast_table.xlsx')

### Step 2: HTTK deduce Css_Plasma with in vitro ####

#Steady State Plasma Concentration Css 
#calc_mc_css(): predict Css resulting from an ongoing 1 mg/kg/day exposure.

for (this.id in unique(toxcast.table$DTXSID)) {
  if (this.id %in% get_cheminfo(info = "dtxsid", suppress.messages = TRUE)) {
    set.seed(12345)
    
    css_results <- calc_mc_css(dtxsid = this.id,
                               output.units = "uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # httk 默认采样数
                               suppress.messages = TRUE)
    
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_5th"]  <- css_results[3]
    
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

for (this.id in unique(toxcast.table$DTXSID)) {
  
  if (this.id %in% get_cheminfo(info="dtxsid", suppress.messages=TRUE) &
      is.na(toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"])) {
    
    set.seed(12345)
    
    css_results <- calc_mc_css(dtxsid=this.id,
                               output.units="uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  
                               suppress.messages=TRUE)
    
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_5th"] <- css_results[3]
    
    toxcast.table[toxcast.table$DTXSID==this.id, "Css.Type"] <- "in silico"
  }
}

#Take another look at our table – many of the gaps are now filled in.

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th',"Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")

sum(!is.na(toxcast.table$Css_50th))

## Output
HBM_Css_Plasma = toxcast.table
HBM_Css_Plasma <- HBM_Css_Plasma %>%
  select(DTXSID, Css_95th, Css_50th, Css_5th, Css.Type) %>%
  dplyr::rename(Css_plasma_95th = Css_95th,
                Css_plasma_50th = Css_50th,
                Css_plasma_5th = Css_5th,
                Css.Type_Plasma = Css.Type)

### Deduce Blood Css ####

### Step 2.1: HTTK deduce Css_Blood with in vitro ####

#Steady State Plasma Concentration Css 
#calc_mc_css(): predict Css resulting from an ongoing 1 mg/kg/day exposure.

for (this.id in unique(toxcast.table$DTXSID)) {
  if (this.id %in% get_cheminfo(info = "dtxsid", suppress.messages = TRUE)) {
    set.seed(12345)
    
    css_results <- calc_mc_css(dtxsid = this.id,
                               output.units = "uM",
                               concentration = 'blood',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  
                               suppress.messages = TRUE)
    
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_5th"]  <- css_results[3]
    
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

for (this.id in unique(toxcast.table$DTXSID)) {
  
  if (this.id %in% get_cheminfo(info="dtxsid", suppress.messages=TRUE) &
      is.na(toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"])) {
    
    set.seed(12345)
    
    css_results <- calc_mc_css(dtxsid=this.id,
                               output.units="uM",
                               concentration = 'blood',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  
                               suppress.messages=TRUE)
    
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_5th"] <- css_results[3]
    
    toxcast.table[toxcast.table$DTXSID==this.id, "Css.Type"] <- "in silico"
  }
}

#Take another look at our table – many of the gaps are now filled in.

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th',"Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")

sum(!is.na(toxcast.table$Css_50th))

## Output
HBM_Css_Blood = toxcast.table
HBM_Css_Blood <- HBM_Css_Blood %>%
  select(DTXSID, Css_95th, Css_50th, Css_5th, Css.Type) %>%
  dplyr::rename(Css_blood_95th = Css_95th,
                Css_blood_50th = Css_50th,
                Css_blood_5th = Css_5th,
                Css.Type_Blood = Css.Type)

##### Merge plasma and blood Css ####
HBM_Css <- left_join(HBM_Css_Blood, HBM_Css_Plasma, by = 'DTXSID')

### Extract predicted internal blood from HExpPredict ####
load("/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/Internal_Prediction.RData")
Predict_Blood <- Internal_Prediction %>%
  select(DTXSID, P5, P50, P95) %>%
  dplyr::rename(CB_P5_uM = P5,
                CB_P50_uM = P50,
                CB_P95_uM = P95) 

HBM_Predict_Blood <- left_join(Che_list, Predict_Blood, by = 'DTXSID')

### Extract predicted external exposure from SEEM3 ####
load('/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/External_Prediction.RData')
Predict_External <- External_Prediction %>%
  select(DTXSID, P5, P50, P95, Concentration_Unit) %>%
  dplyr::rename(Ex_P5 = P5,
                Ex_P50 = P50,
                Ex_P95 = P95,
                Ex_Unit = Concentration_Unit)

HBM_Predict <- left_join(HBM_Predict_Blood, Predict_External, 
                         by = 'DTXSID')

##### Merge ####
Pred_In_Ex <- left_join(HBM_Predict, HBM_Css_Blood, by = 'DTXSID')
Pred_In_Ex_1 <- left_join(Pred_In_Ex, HBM_Css_Plsama, by = 'DTXSID')

library(writexl)
write_xlsx(Pred_In_Ex_1, '5.2-Overlap_HBM_Prediction.xlsx')


####### 6.Plot_HBM_Region_P5.P50.P95--------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

## Load Overlapping compounds HBM_Conc data
library(readxl)
HBM_Conc <- read_excel('2.2-Overlap_HBM_Conc_Unit_Conversion.xlsx',
                       sheet = 1)
unique(HBM_Conc$Unit)

## Add Region_new variable
HBM <- HBM_Conc %>%
  mutate(Region_new = case_when(
    Region == "Northern EU" | Region == "Western EU" | 
      Region == "Southern EU" | Region == "Eastern EU" ~ 'EU',
    Region == "California" | Region ==  "Washington" ~ 'US',
    Country == "United States" ~ 'US',
    Country == "Canada" ~ 'CA')
  )

unique(HBM$Region_new) 
unique(HBM$Sample_INF) 
unique(HBM$Unit)

## Create Unit_new variable: Sample_Unit ##
HBM <- HBM %>%
  mutate(Conc_new = Concentration, 
         Unit_new = paste(Sample_INF, Unit, sep = '_'))

unique(HBM$Unit_new) 

#### Calculate statistic vector p5, p50, p95 for concentration ####
# Updated percentiles function to handle NA values
percentiles <- function(x) {
  quantile(x, probs = c(0.05, 0.50, 0.95), na.rm = TRUE)
}

# Updated summarise code
HBM_1 <- HBM %>%
  group_by(InChIKey, Region_new, Unit_new, Sample_INF, Sample, Sample_Type,
           Chemical_Name, Chemical_Group, DTXSID) %>%
  summarise(
    Conc_P5 = percentiles(Conc_new)[1],
    Conc_P50 = percentiles(Conc_new)[2],
    Conc_P95 = percentiles(Conc_new)[3],
    .groups = "drop"
  ) %>%
  mutate(Region_new = ifelse(Region_new == 'US', 'American', Region_new),
         Region_new = ifelse(Region_new == "CA", 'Canadian', Region_new),
         Region_new = ifelse(Region_new == 'EU', 'European', Region_new))

unique(HBM_1$Unit_new)
unique(HBM_1$InChIKey) #76 compounds

## Filter, Output, and Save data 
# Remove records: Conc_P50 = NA
HBM_2 <- HBM_1 %>%
  filter(!is.na(Conc_P50))

unique(HBM_2$InChIKey) #76 compounds

# Output
library(writexl)
write_xlsx(HBM_2, '6-Overlap_HBM_Conc_P5.P50.P95_N=76.xlsx')


##### Before draw plot: Revise Chemical_Group values and Chemical_Name ####
unique(HBM_2$Chemical_Group)
unique(HBM_2$Chemical_Name)

HBM_3 <- HBM_2 %>%
  # Revise Chemical_Group values
  mutate(Chemical_Group_new = case_when(
    Chemical_Group == "polycyclic aromatic hydrocarbons (PAHs)" ~ "PAHs",
    Chemical_Group == "metals and trace elements" ~ "Metals",
    Chemical_Group == "pesticides: organochlorines"  ~ "Pesticides",
    Chemical_Group == "plasticizers: phthalates" ~ "Phthalates",
    Chemical_Group == "flame retardants" ~ "FRs",
    Chemical_Group == "chlorophenols" ~ "CPs",
    Chemical_Group == "personal care and consumer product chemicals"  ~ "PCPs",
    Chemical_Group == "per and polyfluoroalkyl substances (PFAS)"  ~ "PFAS",
    Chemical_Group == "pesticides: organophosphates" ~ "Pesticides",
    Chemical_Group == "pesticides: pyrethroids" ~ "Pesticides",
    Chemical_Group == "polychlorinated biphenyls (PCBs)" ~ "PCBs",
    Chemical_Group == "nicotine"~ "Nicotine")
  ) %>%
  
  mutate(Chemical_Name = case_when(
    Chemical_Name == "3-(2,2-Dichloroethenyl)-2,2-dimethylcyclopropanecarboxylic acid" ~
      "Permethric acid", 
    TRUE ~ Chemical_Name)
  ) %>%
  
  mutate(Unit_new = case_when(
    Unit_new == "Urine_µg/g" ~ "Urine_µg/g creatinine" ,
    Unit_new == "Breast milk_µg/g breast milk" ~ "Breast milk_µg/g",
    TRUE ~ Unit_new)
  ) 

##### Create plot #####
library(ggplot2)
library(tidyr)
library(dplyr)

HBM_3 <- HBM_3 %>%
  arrange(Chemical_Group_new, Conc_P50) %>%
  mutate(Chemical_Name = factor(Chemical_Name, levels = unique(Chemical_Name)))

unique(HBM_3$Unit_new)
# Define color mapping based on Unit_new
colors <- c("Urine_µmol/L" = '#008000', "Urine_µg/g creatinine" = "#808000", 
            "Urine_µg/L, normalized for SG" = '#77BFC7',
            "Urine_µmol/g creatinine" = '#99C68E', "Urine_µg/L" = '#01F9C6',
            "Urine_µmol/L, normalized for SG" = "#B3D9D9",
            "Urine_µg As/L" = '#08A04B',
            "Urine_µg As/g creatinine" = "#00827F",
            
            "Blood_ng/g serum" = 'red',
            "Blood_µg/L" = "#E77471", 
            "Blood_ng/g lipid" = '#F778A1', 
            
            "Cord blood_µg/L"  = 'orange', "Cord blood_ng/g lipid"  = '#F5E216',
            
            "Breast milk_µg/L" = '#B048B5', "Breast milk_µg/g" = '#B666D2',
            
            "Hair_µg/g" = "black",
            
            "Placenta_µg/g" = 'blue'
)

# Define shape mapping based on Region_new
shapes <- c("European" = 1, "American" = 2, "Canadian" = 4)

### Plot the graph ####
ggplot() +
  # Add median points for HBM_Conc
  geom_point(data = HBM_3, 
             aes(x = Conc_P50, y = reorder(Chemical_Name, Conc_P50),
                 color = Unit_new, shape = Region_new), 
             size = 2.0, alpha = 0.7) +
  
  # Add error bars for 5th to 95th percentiles
  geom_errorbar(data = HBM_3,
                aes(xmin = Conc_P5, xmax = Conc_P95, 
                    y = reorder(Chemical_Name, Conc_P50), color = Unit_new),
                width = 0.1) +
  
  # Set x-axis to log scale
  scale_x_log10(limits = c(1e-5, 1e3),
                breaks = 10^seq(-5, 3, by = 1),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  
  # Add labels
  labs(x = "Concentration range", y = "") +
  
  # Customize theme
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7, hjust = 1),
    axis.text.x = element_text(angle = 0, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.1, "cm"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 5)  # Adjust grouping labels
  ) +
  
  # Custom color legend for Unit_new
  scale_color_manual(name = "Sample Type", values = colors) +
  
  # Custom shape legend for Region_new
  scale_shape_manual(name = "Population", values = shapes) +
  
  # Add grouping
  facet_grid(Chemical_Group_new ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.text.y = element_text(size = 7))  # Adjusts the facet label size







####### 7.Plot_HBM_Pred_Expo_AC50-----------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

library(readxl)
library(dplyr)

## Load HBM_Cytotoxicity AC50 data

Cyto <- read_excel('4.3-Overlap_HBM_Cyto_AC50_Assay.xlsx',
                   sheet = 1)

library(tidyr)

Cyto_1 <- Cyto %>%
  select(-spid, -ac50_2731, -ac50_3070, -ac50_3098, -ac50_1971, 
         -ac50_3163, -ac50_2031) %>%
  rename(DTXSID = dtxsid) %>%
  pivot_longer(cols = starts_with("ac50_"), 
               names_to = "ac50_type", 
               values_to = "ac50_value",
               values_drop_na = TRUE) 

unique(Cyto_1$DTXSID) #45 compounds

## Load HBM_Prediction Data
HBM_Prediction <- read_excel('5.2-Overlap_HBM_Prediction.xlsx',
                             sheet = 1) %>%
  select(DTXSID, Css_blood_5th, Css_blood_50th, Css_blood_95th, Css.Type_Blood,
         CB_P5_uM, CB_P50_uM, CB_P95_uM,
         Ex_P5, Ex_P50, Ex_P95)

## Load HBM_Expo Data ####
HBM_Conc <- read_excel('6-Overlap_HBM_Conc_P5.P50.P95_N=76.xlsx',
                       sheet = 1) %>%
  mutate(Unit = sub(".*_", "", Unit_new))

# Input compound AVERAGE_MASS
Che_List <- read_excel('1-Overlap_EU_US_CA_N=76.xlsx',
                       sheet = 1) %>%
  select(DTXSID, AVERAGE_MASS)

HBM_Conc <- left_join(HBM_Conc, Che_List, by = 'DTXSID')

# Create Sample_Detail
unique(HBM_Conc$Sample_Type)
unique(HBM_Conc$Sample)

HBM_Conc_1 <- HBM_Conc %>%
  mutate(
    Sample_Detail = case_when(
      Sample_INF == "Urine" ~ "Urine",
      
      Sample == "Human (Hair)" ~ "Hair",
      
      Sample == "Placenta tissue (ng/g)" ~ "Placenta",
      Sample == "Breast Milk" ~ "Breast milk",
      
      Sample == "serum" ~ "Blood serum",
      Sample == "Serum" ~ "Blood serum",
      Sample == "plasma" ~ "Blood plasma",
      Sample == "pooled serum" ~ "Blood serum",
      Sample == "blood" ~ "Blood",
      
      Sample_Type == "Blood Whole Blood" ~ "Blood",
      Sample_Type == "Whole Blood" ~ "Blood",
      Sample_Type == "Blood Serum" ~ "Blood serum",
      Sample_Type == "Serum" ~ "Blood serum",
      Sample_Type == "Blood Plasma" ~ "Blood plasma",
      
      
      Sample_Type == "Cord Blood Whole Blood" ~ "Cord blood",
      Sample_Type == "cord whole blood" ~ "Cord blood",
      Sample_Type == "cord plasma fat" ~ "Cord blood plasma fat",
      Sample_Type == "cord plasma" ~ "Cord blood plasma",
      Sample_Type == "Cord Blood Plasma" ~ "Cord blood plasma",
      Sample_Type == "cord plasma fat" ~ "Cord blood plasma",
      Sample_Type == "Cord Blood Serum" ~ "Cord blood serum",
      
      TRUE ~ Sample  # Keep Sample's original value by default
    ))

unique(HBM_Conc_1$Sample_Detail)

unique(HBM_Conc$Unit)
# Convert unit "µg/L" to "µmol/L"
HBM_Conc_2 <- HBM_Conc_1 %>%
  mutate(
    Conc_P5_uM = ifelse(Unit == "µg/L", Conc_P5 / AVERAGE_MASS, NA_real_),
    Conc_P50_uM = ifelse(Unit == "µg/L", Conc_P50 / AVERAGE_MASS, NA_real_),
    Conc_P95_uM = ifelse(Unit == "µg/L", Conc_P95 / AVERAGE_MASS, NA_real_),
    Unit_uM = ifelse(Unit == "µg/L", "µmol/L", NA_character_)
  )

#### HBM_Conc: The function calculates the 5th, 50th, and 95th percentiles of a given vector ####
# Updated percentiles function to handle NA values
percentiles <- function(x) {
  quantile(x, probs = c(0.05, 0.50, 0.95), na.rm = TRUE)
}

# Updated summarise code
HBM_Conc_3 <- HBM_Conc_2 %>%
  select(-Conc_P5, -Conc_P50, -Conc_P95, -Unit, -Unit_new,
         -Conc_P5_uM, -Conc_P95_uM) %>%
  group_by(InChIKey, Sample_Detail, DTXSID,
           Chemical_Name, Chemical_Group) %>%
  summarise(
    P5_Conc_P50_uM = percentiles(Conc_P50_uM)[1],
    P50_Conc_P50_uM = percentiles(Conc_P50_uM)[2],
    P95_Conc_P50_uM = percentiles(Conc_P50_uM)[3],
    .groups = "drop") 

unique(HBM_Conc_3$InChIKey) #76 compounds
unique(HBM_Conc_3$DTXSID)

## Save
library(writexl)
write_xlsx(HBM_Conc_3, '7-Overlap_HBM_P50_uM_P5.P50.P95.xlsx')

##### Merge HBM_Conc, HBM_HBM_Prediction, Cyto_AC50 ####
df <- left_join(HBM_Conc_3, HBM_Prediction, by = 'DTXSID')
df_1 <- left_join(df, Cyto_1, by = 'DTXSID')

unique(df_1$InChIKey) #76 compounds

HBM_Pred_Expo_CytoAC50 <- df_1
save(HBM_Pred_Expo_CytoAC50, file = '7-HBM_Pred_Expo_CytoAC50_N=76.RData')

#### Plot: Cyto_AC50, CB, Css, Expo, SEEM3 ####

load('7-HBM_Pred_Expo_CytoAC50_N=76.RData')
df <- HBM_Pred_Expo_CytoAC50

names(df)
unique(df$Sample_Detail)
unique(df$Chemical_Group)

## Revise Chemical_Group values and Chemical_Name
df_1 <- df %>%
  # Revise Chemical_Group values
  mutate(Chemical_Group_new = case_when(
    Chemical_Group == "polycyclic aromatic hydrocarbons (PAHs)" ~ "PAHs",
    Chemical_Group == "metals and trace elements" ~ "Metals",
    Chemical_Group == "pesticides: organochlorines"  ~ "Pesticides",
    Chemical_Group == "plasticizers: phthalates" ~ "Phthalates",
    Chemical_Group == "flame retardants" ~ "FRs",
    Chemical_Group == "chlorophenols" ~ "CPs",
    Chemical_Group == "personal care and consumer product chemicals"  ~ "PCPs",
    Chemical_Group == "per and polyfluoroalkyl substances (PFAS)"  ~ "PFAS",
    Chemical_Group == "pesticides: organophosphates" ~ "Pesticides",
    Chemical_Group == "pesticides: pyrethroids" ~ "Pesticides",
    Chemical_Group == "polychlorinated biphenyls (PCBs)" ~ "PCBs",
    Chemical_Group == "nicotine"~ "Nicotine")
  ) %>%
  
  mutate(Chemical_Name = case_when(
    Chemical_Name == "3-(2,2-Dichloroethenyl)-2,2-dimethylcyclopropanecarboxylic acid" ~
      "Permethric acid", 
    TRUE ~ Chemical_Name)
  ) 

unique(df_1$Chemical_Group_new)

library(ggplot2)
library(patchwork)
# Ensure that all compounds are visible on the y-axis and grouped by Chemical_Group
df_1 <- df_1 %>%
  arrange(Chemical_Group, Ex_P50) %>%
  mutate(Chemical_Name = factor(Chemical_Name, levels = unique(Chemical_Name)))

## Plot 1: CB, Css, Exposome (Conc_uM) ####
plot1 <- ggplot(df_1, aes(y = reorder(Chemical_Name, Ex_P50))) +  
  
  # Add AC50 points
  geom_point(aes(x = ac50_value, color = "Cytotoxicity_AC50"), 
             size = 1, shape = 3, alpha = 0.7) +
  
  # Add CB error bars to represent the range from 5th to 95th percentiles
  geom_errorbar(aes(xmin = CB_P5_uM, 
                    xmax = CB_P95_uM, color = "CB_95%CI"), 
                width = 0.1) +
  # 
  geom_point(aes(x = CB_P50_uM, color = "CB_Median"), 
             size = 1.5, shape = 16) +
  
  # Add Css_blood error bars to represent the 5th to 95th percentile range
  geom_errorbar(aes(xmin = Css_blood_5th, 
                    xmax = Css_blood_95th, color = "Css_blood_95%CI"), 
                width = 0.1) +
  
  # Add Css_blood median point
  geom_point(aes(x = Css_blood_50th, color = "Css_blood_Median"), 
             size = 1.5, shape = 16, alpha = 0.7) +
  
  # Add Conc error bars
  geom_errorbar(aes(xmin = P5_Conc_P50_uM, 
                    xmax = P95_Conc_P50_uM, color = "Measured_Conc_95%CI"), 
                width = 0.1) +
  
  # Add Conc median point and group different shapes according to Conc_Sample
  geom_point(aes(x = P50_Conc_P50_uM, color = "Measured_Conc_Median", 
                 shape = Sample_Detail), size = 1.5, alpha = 0.7) +
  
  # Set the x-axis to logarithmic scale
  scale_x_log10(limits = c(1e-7, 1e5),
                breaks = 10^seq(-7, 5, by = 1),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  
  # Add tags
  labs(x = "Concentration (µM)", color = "Legend", shape = "Sample Type") +
  
  # Custom colors
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
                                "Breast milk" = 4, 
                                "Cord blood" = 2,
                                "Cord blood plasma" = 6,
                                "Cord blood serum" = 7,
                                "Hair" = 3,
                                "Blood serum" = 5,
                                "Placenta" = 9)
  ) + 
  
  # Remove the y-axis label and ticks
  theme_minimal() +
  theme(axis.title.y = element_blank(),  # 
        axis.text.y = element_blank(),   # 
        axis.ticks.y = element_blank(),  #
        axis.text.x = element_text(size = 7)) +
  
  # Add grouping
  facet_grid(Chemical_Group_new ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.text.y = element_text(size = 9))  # Adjusts the facet label size


print(plot1)


## Plot 2: External Prediction SEEM3 ####
plot2 <- ggplot(df_1, aes(y = reorder(Chemical_Name, Ex_P50))) +  
  
  # Added SEEM3 error bars to represent the range of 5th to 95th percentiles
  geom_errorbar(aes(xmin = Ex_P5, 
                    xmax = Ex_P95, color = "Exposure_95%CI"), 
                width = 0.1) +
  
  # Added SEEM3 median point
  geom_point(aes(x = Ex_P50, color = "Exposure_Median"), 
             size = 1.5, shape = 0) +  
  
  # Set the x-axis to logarithmic scale
  scale_x_log10(limits = c(1e-21, 1e2),
                breaks = 10^seq(-21, 2, by = 2),
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
        axis.text.x = element_text(size = 7.25)) +
  
  # Add grouping
  facet_grid(Chemical_Group_new ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.text.y = element_text(size = 9))  # Adjusts the facet label size


print(plot2)

## Merge two charts, align them horizontally, and share a common y-axis ####
combined_plot <- plot2 + plot1 + 
  plot_layout(ncol = 2, guides = "collect") & theme(legend.position = 'right')

# Output the merged chart
combined_plot

















####### 8.HBM_LD50_ChemPro_CTV_Cytotoxicity-------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

##### 8.1. ChemPro and LD50 data for HBM 76 compounds ####
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

##### 8.2. Conditional Toxicity Value (CTV) #####
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

##### 8.3. CompTox Cytotoxicity Data #####

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

library(tidyr)
library(dplyr)

HBM_Cyto_AC50_unique <- HBM_Cyto_AC50 %>%
  distinct(dtxsid, aeid, ac50, .keep_all = TRUE) #3675 records

HBM_Cyto_AC50_unique <- HBM_Cyto_AC50_unique %>%
  group_by(dtxsid, aeid, spid) %>%
  mutate(record_id = row_number()) %>%
  ungroup()

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

HBM_Cyto_HIT_unique <- HBM_Cyto_HIT %>%
  distinct(dtxsid, aeid, HIT_CALL, .keep_all = TRUE) #2835 records

HBM_Cyto_HIT_unique <- HBM_Cyto_HIT_unique %>%
  group_by(dtxsid, aeid, spid) %>%
  mutate(record_id = row_number()) %>%
  ungroup()

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

HBM_Cyto_HIT_Assay_numeric <- HBM_Cyto_HIT_Assay %>%
  mutate(across(starts_with("HIT"), ~ ifelse(. == "Active", 1, 0)))

str(HBM_Cyto_HIT_Assay_numeric)
print(HBM_Cyto_HIT_Assay_numeric)


HBM_Cyto_HIT_Assay_percentage <- HBM_Cyto_HIT_Assay_numeric %>%
  group_by(dtxsid) %>%
  summarise(
    total_hits = sum(!is.na(across(starts_with("HIT_")))), 
    active_count = sum(across(starts_with("HIT_"), ~ sum(. == 1, na.rm = TRUE))), 
    ratio = (active_count / total_hits) * 100, 
    .groups = 'drop' 
  )

print(HBM_Cyto_HIT_Assay_percentage) # 45 records

HBM_Cyto_HIT_Assay_Ratio <- HBM_Cyto_HIT_Assay %>%
  left_join(HBM_Cyto_HIT_Assay_percentage, by = "dtxsid") %>%
  distinct(dtxsid, .keep_all = TRUE) #45 records

print(HBM_Cyto_HIT_Assay_Ratio)

save(HBM_Cyto_HIT_Assay_Ratio, file = '8.6-HBM_Cyto_%ActiveAssay_N=45.RData')

#### Output Cytoxicity maxMeanConc #####
# Extract AEID, maxMeanConc, DTXSID
HBM_Cyto_Max <- subset(HBM_Cyto_2,
                       select = c(dtxsid, aeid, spid, maxMeanConc)) #3944 records

library(tidyr)
library(dplyr)

HBM_Cyto_Max_unique <- HBM_Cyto_Max %>%
  distinct(dtxsid, aeid, maxMeanConc, .keep_all = TRUE) #3709 records

HBM_Cyto_Max_unique <- HBM_Cyto_Max_unique %>%
  group_by(dtxsid, aeid, spid) %>%
  mutate(record_id = row_number()) %>%
  ungroup()

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



####### 9.ToxPi_HBM_Data_Preparation--------------------------------------------

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












###### 10.ToxPi_Plot_HBM--------------------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

#### Toxpi_score: HBM_Overlapping compounds
load('9-HBM_Tox_N=76.RData') #76 compounds

library(dplyr)
## Remove all metals
HBM_Tox <- HBM_Tox %>%
  filter(Chemical_Group != 'metals and trace elements') %>%
  rename(LD50 = `ORAL_RAT_LD50_MOL/KG_TEST_PRED`,
         LogKoa = `LogKoa: Octanol-Air`,
         LogKow = `LogKow: Octanol-Water`) %>%
  
  mutate(PREFERRED_NAME = case_when(
    PREFERRED_NAME == "3-(2,2-Dichloroethenyl)-2,2-dimethylcyclopropanecarboxylic acid" ~
      "Permethric acid", 
    TRUE ~ PREFERRED_NAME)
  ) #62 compounds

names(HBM_Tox)
# Find the location of the hit variable
hit_index <- which(names(HBM_Tox) == "HIT_Ratio")
# Extract all column names after the hit variable
ac50 <- names(HBM_Tox)[(hit_index + 1):ncol(HBM_Tox)]

#### Toxpi model setting ####
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

# Automatically extract all column names starting with 'ac50_'
ac50_columns <- grep("^ac50_", names(HBM_Tox), value = TRUE)

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
  Slice12 = TxpSlice('LD50')
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

final.trans <- list(
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
  }
)

f.model <- TxpModel(txpSlices = f.slices,
                    txpWeights = c (1,1,1,1,1,1,1,1,1,1,1,1),
                    txpTransFuncs = final.trans)

## Calculate ToxPi scores

f.results <- txpCalculateScores(model = f.model,
                                input = HBM_Tox,
                                id.var = 'PREFERRED_NAME')
# ToxPi scores
toxpi_scores <- txpSliceScores(f.results)
print(toxpi_scores)


#### Pie Plotting ####

## Setting colors
library(RColorBrewer)
# Using RColorBrewer palettes
colors <- brewer.pal(n = 12, name = "Paired") # Color palette of up to 12 colors

## Draw plot
library(grid)
plot(f.results,
     showScore = TRUE,
     fills = colors # fill colors
)

#### Edit Pie Plotting ####
# List all grobs in the current drawing
grid.ls()
# Find all grobs containing 'label' and change the font size
all_grobs <- grid.ls(print = FALSE)


## Modify the font size of all labels
label_grobs <- grep("label", all_grobs$name, value = TRUE)

for (label in label_grobs) {
  grid.edit(label, gp = gpar(fontsize = 6))
}

## Edit 'Toxpi Scores' characteristics
for (s in sprintf("pie-%d-radSum", 1:63)) {
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

HBM_Tox_1 <- HBM_Tox %>%
  mutate(id = row_number()) %>%
  select(id, DTXSID, PREFERRED_NAME)

# Transform matrix to dataframe
toxpi_scores <- as.data.frame(toxpi_scores) 

# Add toxpi_score variable which is the sum of all slices
toxpi_scores <- toxpi_scores %>%
  mutate(id = row_number()) %>%
  rowwise() %>%
  mutate(toxpi_score = sum(c_across(Slice1:Slice12))) %>%
  ungroup()

# Merged chemcial info and slice, toxpi scores
HBM_Toxpi_Slices <- HBM_Tox_1 %>%
  full_join(toxpi_scores, by = 'id')

# Output ToxpiScore
library(writexl)
write_xlsx(HBM_Toxpi_Slices,'10-HBM_Tox_Slice_N=63.xlsx')

#### Upgrade Rank plot ####
library(ggplot2)

RankPlot <- ggplot(HBM_Toxpi_Slices, 
                   aes(x = toxpi_score, y = reorder(PREFERRED_NAME, toxpi_score))) +
  geom_point() +
  scale_x_continuous(limits = c(0.30, 0.65), breaks = seq(0.30, 0.65, 0.05)) +
  labs(x = "Tox Score", y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(
      hjust = 1,  # 
      size = 6.0    #
    )
  )

show(RankPlot)

# Output plot
ggsave("10.2-RankPlot_HBM_Tox.pdf", plot = RankPlot, 
       width = 6, height = 8, units = "in", dpi = 600)


### Hierarchical Clustering ####
f.hc <- hclust(dist(txpSliceScores(f.results)))
plot(f.hc, hang = -1, labels = txpIDs(f.results), xlab="", sub="")








####### 11.ExpTox_HBM_BEQ_BER_Calculation---------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

library(readxl)
library(dplyr)

#### API-CompTox-ctxR package, extract Target assay data ####

library(ctxR)
### API key will be stored as the variable my_key
my_key <- '5e944183-665f-4e1b-9ff0-c7e27b8b09e4'

### Retrieve all assays
assays <- get_all_assays(API_key = my_key) #1499 assays

library(writexl)
write_xlsx(assays, 'CompTox_assays.xlsx')

# Load Target Assays
TargetAssay <- read_excel('TargetedAssays.xlsx',
                          sheet = "HBM")

# Extract AEID for Target Assays and convert to Values
TA_AEID <- TargetAssay %>%
  select(aeid) %>%
  distinct() %>%
  pull(aeid) #12 AEIDs

# Extract DTXSID from HBM 
HBM <- read_excel('1-Overlap_EU_US_CA_N=76.xlsx',
                  sheet = 1)
DTXSID <- HBM %>%
  select(DTXSID) %>%
  distinct() %>%
  pull(DTXSID) #76 DTXSID

### Retrieve bioactivity data from NR_AEID batch
TA_bio_HBM <- get_bioactivity_details_batch(AEID = c(785,786,761,764,762,1816,
                                                     802,801,1198,1128,804,805), #12 Target Assays
                                            API_key = my_key,
                                            DTXSID = DTXSID #HBM compounds
)
save(TA_bio_HBM, file = "11.1-HBM_TargetAssay_Bioactivity.RData")

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
load("11.1-HBM_TargetAssay_Bioactivity.RData")
TA_bio_HBM_filtered <- map(TA_bio_HBM, ~ {
  if (is_valid_df(.x)) .x else NULL
})

# Remove NULL elements from the list
TA_bio_HBM_filtered <- compact(TA_bio_HBM_filtered)

TA_bio_HBM_merge <- TA_bio_HBM_filtered %>%
  map_df(~ .x) #33480 records

unique(TA_bio_HBM_merge$dtxsid) #48 compounds

unique(TA_bio_HBM_merge$aeid)

# Save data
save(TA_bio_HBM_merge, file = '11.2-HBM_TA_Bioactivity_N=48.RData')

unique(TA_bio_HBM_merge$dtxsid) #48 compounds

##### Calculate BEQ ####
## Select target aeid 
load("11.2-HBM_TA_Bioactivity_N=48.RData")
HBM_TA_bio <- TA_bio_HBM_merge %>%
  semi_join(TargetAssay, by = 'aeid') #619 records

unique(HBM_TA_bio$aeid) #12 aeid
unique(HBM_TA_bio$dtxsid) #36 compounds

HBM_TA_bio_1 <- HBM_TA_bio %>%
  select(dtxsid, aeid, ac50) %>%
  group_by(dtxsid, aeid) %>%
  filter(ac50 ==  min(ac50)) %>%
  sample_n(1) %>%
  ungroup()

## calculate ac50_ref for assay: select min ac50 -------------------------------

ac50_ref_df <- HBM_TA_bio_1 %>%
  group_by(aeid) %>%
  dplyr::summarise(ac50_ref_min = min(ac50, na.rm = TRUE)) %>%
  ungroup()

## Load HBM_CB data ---------------------------------------------
load('7-HBM_Pred_Expo_CytoAC50_N=76.RData')
unique(HBM_Pred_Expo_CytoAC50$DTXSID)

HBM_CB <- HBM_Pred_Expo_CytoAC50 %>%
  select(DTXSID, CB_P50_uM) %>%
  rename(dtxsid = DTXSID) %>%
  distinct()

# merge CB and ac50_ref
df_base <- HBM_TA_bio_1 %>%
  left_join(HBM_CB,  by = "dtxsid") %>%
  left_join(ac50_ref_df, by = "aeid") %>%
  distinct()


## Calculate CB/AC50, BEQ, BEQ% -----------------------------------
df_beq <- df_base %>%
  mutate(CB_AC50 = CB_P50_uM / ac50) %>%
  mutate(BEQ_min = CB_AC50 * ac50_ref_min) %>%
  group_by(aeid) %>%
  mutate(sum_BEQ_min = sum(BEQ_min, na.rm = TRUE),
         BEQ_percent_min = BEQ_min / sum_BEQ_min * 100) %>%
  ungroup()

# Each compound corresponding to multiple target assays
library(tidyr)
df_beq_1 <- df_beq %>%
  select(dtxsid, aeid, BEQ_percent_min) %>%
  pivot_wider(names_from = aeid,
              values_from = BEQ_percent_min) 

HBM_TA_BEQ_CB <- df_beq_1
save(HBM_TA_BEQ_CB, file = '11.3-HBM_TA_BEQ_CB_N=36.RData')

##### Calculate BER #####

## Load HBM Overlapping compounds N=76 ##
library(readxl)
Che_list <- read_excel('1-Overlap_EU_US_CA_N=76.xlsx',
                       sheet = 1)

### HTTK model to calculate Css ###
library(httk)
library(readxl)
library(ggplot2)

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
write_xlsx(toxcast.table, '11-HBM_HTTK_toxcast_table.xlsx')

### Step 2: HTTK deduce Css_Plasma with in vitro ####

#Steady State Plasma Concentration Css 
#calc_mc_css(): predict Css resulting from an ongoing 1 mg/kg/day exposure.

# Assume toxcast.table already has columns "Css_95th", "Css_50th", "Css_5th"
for (this.id in unique(toxcast.table$DTXSID)) {
  # Check if a chemical substance is in the HTTK database
  if (this.id %in% get_cheminfo(info = "dtxsid", suppress.messages = TRUE)) {
    # Set the random number seed to ensure reproducibility of Monte Carlo simulations
    set.seed(12345)
    
    # Calculating CSS at different percentiles
    css_results <- calc_mc_css(dtxsid = this.id,
                               output.units = "uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  
                               suppress.messages = TRUE)
    
    # Store the results in a table
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_5th"]  <- css_results[3]
    
    # Mark CSS Type
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

# Make sure the Css column exists, if not add it
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
    
    # 
    set.seed(12345)
    
    # Calculate 95th, 50th and 5th percentile CSS
    css_results <- calc_mc_css(dtxsid=this.id,
                               output.units="uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  
                               suppress.messages=TRUE)
    
    # Store calculation results in different columns
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_5th"] <- css_results[3]
    
    # Mark CSS type as "in silico"
    toxcast.table[toxcast.table$DTXSID==this.id, "Css.Type"] <- "in silico"
  }
}

#Take another look at our table – many of the gaps are now filled in.

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th',"Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")

sum(!is.na(toxcast.table$Css_50th))

## Output
HBM_Css_Plsama = toxcast.table

## Step 4: Reverse Dosimetry In Vitro-In Vivo Extropolation ####

#Don’t forget the ToxCast concentrations are given on the log-10 scale:
HBM_Css_Plsama$EquivDose <- signif(10^HBM_Css_Plsama$Q10.AC50 / HBM_Css_Plsama$Css_50th,
                                   3)

knitr::kable(HBM_Css_Plsama[1:10,c("Compound","Q10.AC50","Css_50th","EquivDose")], 
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
                                     DTXSID %in% HBM_Css_Plsama$DTXSID))

# Merge database
HBM_Css_Plsama <- left_join(HBM_Css_Plsama, example.seem, by = 'DTXSID')

#Calculate a bioactivity:exposure ratio (BER)
HBM_Css_Plsama$BER <- signif(HBM_Css_Plsama$EquivDose/HBM_Css_Plsama$P95,3)

HBM_BER <- HBM_Css_Plsama %>%
  select(InChIKey, DTXSID, BER) %>%
  filter(!is.na(BER))

unique(HBM_BER$InChIKey)
unique(HBM_BER$DTXSID)

save(HBM_BER, file = '11.4-HBM_HTTK_BER_N=20.RData')




####### 12.ExpTox_Plot_HBM------------------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

### BEQ_CB TargetAssay ####
load('11.3-HBM_TA_BEQ_CB_N=36.RData')
BEQ <- HBM_TA_BEQ_CB %>%
  rename(DTXSID = dtxsid) %>%
  rename_with(~ paste0("BEQ_", .), .cols = -1) %>%
  # remove compounds with many NA values
  filter(!(DTXSID %in% c('DTXSID7020508', 'DTXSID9020376'))) #34 compounds

### BER ####
load('11.4-HBM_HTTK_BER_N=20.RData')
BER <- HBM_BER %>%
  select(DTXSID, BER) #20 compounds

### RI calculation using Urine Conc ####
library(readxl)
Toxpi <- read_excel('10-HBM_Tox_Slice_N=63.xlsx',
                    sheet = 1)

HBM <- read_excel('7-Overlap_HBM_P50_uM_P5.P50.P95.xlsx',
                  sheet = 1) %>%
  filter(Sample_Detail == 'Urine') %>%
  filter(Chemical_Group != 'metals and trace elements') %>%
  distinct() 

load('7-HBM_Pred_Expo_CytoAC50_N=76.RData')

unique(HBM$InChIKey) #46
unique(HBM$DTXSID)

HBM_Urine <- HBM %>%
  select(DTXSID, P50_Conc_P50_uM) %>%
  mutate(log_P50 = log(P50_Conc_P50_uM)) # log-transformed concentration

df <- left_join(Toxpi, HBM_Urine, by = 'DTXSID')

## Normalized log_P50

df$Conc_Urine_normalized <- (df$log_P50 - min(df$log_P50, na.rm = TRUE)) / 
  (max(df$log_P50, na.rm = TRUE) - min(df$log_P50, na.rm = TRUE))

## Calculate RI
df <- df %>%
  mutate(RI_Urine = Conc_Urine_normalized * toxpi_score)

#### Merge RI 和 BEQ, BER ####
df_1 <- df %>%
  left_join(BER, by = 'DTXSID') %>%
  left_join(BEQ, by = 'DTXSID') %>%
  
  mutate(PREFERRED_NAME = case_when(
    PREFERRED_NAME == "3-(2,2-Dichloroethenyl)-2,2-dimethylcyclopropanecarboxylic acid" ~
      "Permethric acid", 
    TRUE ~ PREFERRED_NAME)
  )

ExpTox = df_1
save(ExpTox, file = '12-ExpTox_Plot_HBM.RData')

unique(df_1$DTXSID)
unique(df_1$RI_Urine)
unique(df_1$BER)

#### ExpTox calculation ####
library(toxpiR)

# View all functions
lsf.str("package:toxpiR")

## Create slice with transformation functions
# Slice 1: BEQ, # log10(x)-log10(min)
# Slice 2: BER, # log10(x)+log10(max)
# Slice 3: RI, # log10(x)-log10(min)

# Automatically extract all column names starting with 'BEQ_'
BEQ_columns <- grep("^BEQ_", names(df_1), value = TRUE)

f.slices <- TxpSliceList(
  Slice1 = TxpSlice(BEQ_columns),
  Slice2 = TxpSlice('BER'),
  Slice3 = TxpSlice('RI_Urine')
)

## Goal - Create ToxPi model.
# All Slices, weight = 1.

## Method#1: Not consider missing value

final.trans <- TxpTransFuncList(
  f1 = function(x) log10(x),
  f2 = function(x) log10(x),
  f3 = function(x) log10(x))

f.model <- TxpModel(txpSlices = f.slices,
                    txpWeights = c (1,1,1),
                    txpTransFuncs = final.trans)


## Calculate ToxPi scores

f.results <- txpCalculateScores(model = f.model,
                                input = df_1,
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
# List all grobs in the current drawing
grid.ls()
# Find all grobs containing 'label' and change the font size
all_grobs <- grid.ls(print = FALSE)


## Modify the font size of all labels
label_grobs <- grep("label", all_grobs$name, value = TRUE)

for (label in label_grobs) {
  grid.edit(label, gp = gpar(fontsize = 6))
}

## Edit 'Toxpi Scores' characteristics
for (s in sprintf("pie-%d-radSum", 1:63)) {
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
df_1 <- df_1 %>%
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
HBM_ExpTox_Slices <- df_2 %>%
  full_join(toxpi_scores, by = 'id')

library(writexl)
write_xlsx(HBM_ExpTox_Slices, '12-HBM_ExpTox_Slice_N=63.xlsx')

#### Upgrade Rank plot ####
library(ggplot2)

# Remove missing values
HBM_ExpTox_Slices <- HBM_ExpTox_Slices %>%
  filter(!is.na(toxpi_score))

RankPlot <- ggplot(HBM_ExpTox_Slices, 
                   aes(x = toxpi_score, y = reorder(PREFERRED_NAME, toxpi_score))) +
  geom_point() +
  scale_x_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.1)) +
  labs(x = "ExpTox Score", y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(
      hjust = 1,  
      size = 6.0    
    )
  )

show(RankPlot)

# Output plot
ggsave("12.2-RankPlot_HBM_ExpTox.pdf", plot = RankPlot, 
       width = 6, height = 8, units = "in", dpi = 600)



### Hierarchical Clustering ####
f.hc <- hclust(dist(txpSliceScores(f.results)))
plot(f.hc, hang = -1, labels = txpIDs(f.results), xlab="", sub="")







####### 13.Plot_Rank_Change_HBM-------------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

library(readxl)
library(dplyr)

Tox <- read_excel('10-HBM_Tox_Slice_N=63.xlsx',
                  sheet = 1)
Tox_1 <- Tox %>%
  select(PREFERRED_NAME, toxpi_score) %>%
  mutate(rank_tox = rank(-toxpi_score, ties.method = "average")) %>%
  rename(Tox_score = toxpi_score)

ExpTox <- read_excel('12-HBM_ExpTox_Slice_N=63.xlsx',
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
  geom_text_repel(size = 3.5, max.overlaps = Inf) + 
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


####### 14.HBM_TableS4----------------------------------------------------------

setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

library(readxl)
library(dplyr)

df <- read_excel('1-Overlap_EU_US_CA_N=76.xlsx',
                 sheet = 1)

df_1 <- df %>%
  select(InChIKey, DTXSID, PREFERRED_NAME, CASRN, Chemical_Group)

## Load Actual exposure of HBM concentration_uM data
HBM_uM <- read_excel('7-Overlap_HBM_P50_uM_P5.P50.P95.xlsx',
                     sheet = 1)
unique(HBM_uM$Sample_Detail)
HBM_uM_1 <- HBM_uM %>%
  select(InChIKey, Sample_Detail, 
         P5_Conc_P50_uM, P50_Conc_P50_uM, P95_Conc_P50_uM)

## Merge concentration value as P50(P5, P95)
library(stringr)

HBM_uM_2 <- HBM_uM_1 %>%
  mutate(
    Conc_uM = case_when(
      !is.na(P5_Conc_P50_uM) & !is.na(P50_Conc_P50_uM) & !is.na(P95_Conc_P50_uM) ~
        str_c(
          formatC(P50_Conc_P50_uM, format = "e", digits = 3), "(", 
          formatC(P5_Conc_P50_uM, format = "e", digits = 3), ", ", 
          formatC(P95_Conc_P50_uM, format = "e", digits = 3), ")"
        ),
      TRUE ~ NA_character_
    )
  )

## Long form to wide format
library(tidyr)
HBM_uM_2 <- HBM_uM_2[-140, ]
HBM_uM_3 <- HBM_uM_2 %>%
  select(InChIKey, Sample_Detail, Conc_uM) %>%
  pivot_wider(
    id_cols = InChIKey,
    names_from = Sample_Detail,
    values_from = Conc_uM,
    values_fn = list(Conc_uM = ~ paste(unique(.), collapse = ",")),
    values_fill = list(Conc_uM = NA_character_)
  ) 

colnames(HBM_uM_3)
HBM_uM_4 <- HBM_uM_3 %>%
  select(-Hair, -Placenta, -`Cord blood plasma fat`) 

###
df_2 <- left_join(df_1, HBM_uM_4, by = 'InChIKey')

### Load HBM_Compounds Prediction data
HBM_Pred <- read_excel('5.2-Overlap_HBM_Prediction.xlsx',
                       sheet = 1)
HBM_Pred_1 <- HBM_Pred %>%
  select(InChIKey, CB_P5_uM, CB_P50_uM, CB_P95_uM,
         Css_blood_5th, Css_blood_50th, Css_blood_95th,
         Css_plasma_5th, Css_plasma_50th, Css_plasma_95th,
         Ex_P5, Ex_P50, Ex_P95) 

###
df_3 <- left_join(df_2, HBM_Pred_1, by = 'InChIKey')

### Load BEQ% data
load('/Users/ydai/Desktop/Exposome_NGRA/2505-HBM_Analysis/11.3-HBM_TA_BEQ_CB_N=36.RData')
HBM_TA_BEQ_CB_1 <- HBM_TA_BEQ_CB %>%
  rename(DTXSID = dtxsid) %>%
  rename(TOX21_Era_BLA_Agonist_ratio = "785" ,
         TOX21_Era_BLA_Antagonist_ratio = "786",
         TOX21_AR_BLA_Agonist_ratio = "761",
         Tox21_AR_LUC_MDAKB2_Agonist = "764",
         TOX21_AR_BLA_Antagonist_ratio = "762",
         TOX21_AR_LUC_MDAKB2_Antagonist_0.5nM_R1881 = "1816",
         Tox21_PPARg_BLA_Agonist_ratio = "802",
         TOX21_PPARg_BLA_Agonist_ch2 = "801",
         TOX21_PPARg_BLA_Antagonist_ch1 = "1198",
         TOX21_PPARg_BLA_antagonist_viability = "1128",
         TOX21_TR_LUC_GH3_Agonist = "804",
         TOX21_TR_LUC_GH3_Antagonist = "805")

###
df_4 <- left_join(df_3, HBM_TA_BEQ_CB_1, by = 'DTXSID')

### Load BER, RI, Tox
load('12-ExpTox_Plot_HBM.RData')
ExpTox_1 <- ExpTox %>%
  select(DTXSID, toxpi_score, RI_Urine, BER) %>%
  rename(Tox_score = toxpi_score)

###
df_5 <- left_join(df_4, ExpTox_1, by = 'DTXSID')

### Load ExpTox
ExpTox_score <- read_excel('12-HBM_ExpTox_Slice_N=63.xlsx',
                           sheet = 1)
ExpTox_score_1 <- ExpTox_score %>%
  select(DTXSID, toxpi_score) %>%
  rename(ExpTox_score = toxpi_score)

###
df_6 <- left_join(df_5, ExpTox_score_1, by = 'DTXSID')

library(writexl)
write_xlsx(df_6, '14-HBM_Table S4.xlsx')
