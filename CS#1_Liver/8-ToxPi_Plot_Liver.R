
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

#### Toxpi_score: Liver
load('7-Liver_Tox_N=94.RData')

###### Toxpi model conctruction ####
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
# Slice 13: NTP_Group (In vivo evidence), # linear
# Slice 14: IARC, #linear

# 自动提取以 'ac50_' 开头的所有列名
ac50_columns <- grep("^ac50_", names(Liver_Tox), value = TRUE)

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
  Slice12 = TxpSlice('LD50'),
  Slice13 = TxpSlice('NTP_Group'),
  Slice14 = TxpSlice('IARC_Group.y')
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

final.trans <- TxpTransFuncList(
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
  },
  f13 = linear,
  f14 = linear
)


f.model <- TxpModel(txpSlices = f.slices,
                    txpWeights = c (1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                    txpTransFuncs = final.trans)

## Calculate ToxPi scores

f.results <- txpCalculateScores(model = f.model,
                                input = Liver_Tox,
                                id.var = 'PREFERRED_NAME.x')
# ToxPi scores
toxpi_scores <- txpSliceScores(f.results)
print(toxpi_scores)

#### Pie Plotting ####

## Setting colors
library(RColorBrewer)
# Get the 12 colors from the "Paired" palette
colors <- brewer.pal(12, "Paired")
# Add two additional distinct colors
colors <- c(colors, "#FF69B4", "#8A2BE2")

## Draw plot
library(grid)
plot(f.results,
     showScore = TRUE,
     fills = colors # fill colors
)

#### Edit Pie Plotting ####
# 列出当前图形中的所有 grobs
grid.ls()
# 找到所有包含 'label' 的 grobs 并修改字体大小
all_grobs <- grid.ls(print = FALSE)


## 修改所有 label 的字体大小
label_grobs <- grep("label", all_grobs$name, value = TRUE)

for (label in label_grobs) {
  grid.edit(label, gp = gpar(fontsize = 6))
}

## Edit 'Toxpi Scores' characteristics
for (s in sprintf("pie-%d-radSum", 1:94)) {
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

Liver_Tox_1 <- Liver_Tox %>%
  mutate(id = row_number()) %>%
  select(id, DTXSID, PREFERRED_NAME.x)

# Transform matrix to dataframe
toxpi_scores <- as.data.frame(toxpi_scores) 

# Add toxpi_score variable which is the sum of all slices
toxpi_scores <- toxpi_scores %>%
  mutate(id = row_number()) %>%
  rowwise() %>%
  mutate(toxpi_score = sum(c_across(Slice1:Slice14))) %>%
  ungroup()

# Merged chemcial info and slice, toxpi scores
Liver_Toxpi_Slices <- Liver_Tox_1 %>%
  full_join(toxpi_scores, by = 'id') %>%
  rename(PREFERRED_NAME = PREFERRED_NAME.x)

# Output ToxpiScore
library(writexl)
write_xlsx(Liver_Toxpi_Slices,'8-Liver_Tox_Slices_N=94.xlsx')

#### Upgrade Rank plot ####
library(ggplot2)

RankPlot <- ggplot(Liver_Toxpi_Slices, 
                   aes(x = toxpi_score, y = reorder(PREFERRED_NAME, toxpi_score))) +
  geom_point() +
  scale_x_continuous(limits = c(0.2, 0.75), breaks = seq(0.2, 0.75, 0.05)) +
  labs(x = "Tox Score", y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(
      hjust = 1,  # 右对齐
      size = 6.0    # 减小字体大小，可以根据需要调整这个值
    )
  )

show(RankPlot)

# Output plot
ggsave("8.2_RankPlot_Tox_Liver.pdf", plot = RankPlot, 
       width = 6, height = 8, units = "in", dpi = 600)



### Hierarchical Clustering ####
f.hc <- hclust(dist(txpSliceScores(f.results)))
plot(f.hc, hang = -1, labels = txpIDs(f.results), xlab="", sub="")
