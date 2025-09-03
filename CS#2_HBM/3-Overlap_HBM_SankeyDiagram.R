
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
