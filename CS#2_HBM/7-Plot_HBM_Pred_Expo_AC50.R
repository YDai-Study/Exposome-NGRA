
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
  scale_shape_manual(values = c("Urine" = 17,  # 三角形
                                "Blood" = 16,  # 实心圆
                                "Blood plasma" = 1,  # 空心圆
                                "Breast milk" = 4, # 十字形
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














