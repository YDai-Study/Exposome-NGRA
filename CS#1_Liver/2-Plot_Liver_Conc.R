
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')


library(readxl)

Liver <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                    sheet = 1)

# Load Liver_Conc data
Liver_Conc <- read_excel('1.2-Liver_Expo_Conc.xlsx',
                         sheet = 1)
unique(Liver_Conc$InChIKey)
##### 1. Unit conversion ####

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
      TRUE ~ Unit)  # 保留其他 Unit 不变
  )

unique(Liver_Conc_1$Unit)

##### 2. Sample conversion #####
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

