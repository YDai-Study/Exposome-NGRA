
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

# Load All_HBM_Merge data
load('/Users/ydai/Desktop/2503-Exposome_Database_Clean/4-HBM_DataMerge&Clean/All_HBM_Merge_New.RData')
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

##### 1.Overlap Results in EU,US,CA HBM Databases ######
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



##### 2. DonutChart_HBM_ChemGroup_N=76 ####
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




