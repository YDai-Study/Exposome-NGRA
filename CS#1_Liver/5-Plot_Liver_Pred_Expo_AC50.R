
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

## Load Liver_Cytotoxicity AC50 data

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

## Load Liver_Prediction Data
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
  scale_shape_manual(values = c("Urine" = 17,  # 三角形
                                "Blood" = 16,  # 实心圆
                                "Blood plasma" = 1,  # 空心圆
                                "Blood serum" = 5,
                                "Cord blood plasma" = 6,
                                "Plasma" = 1,  # 空心圆
                                "Breast milk" = 4, # 十字形
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











