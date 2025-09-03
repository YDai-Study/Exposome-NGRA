
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

# 确保所有化合物在 y 轴上都能显示，并按 Chemical_Group 分组
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




