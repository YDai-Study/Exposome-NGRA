
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





