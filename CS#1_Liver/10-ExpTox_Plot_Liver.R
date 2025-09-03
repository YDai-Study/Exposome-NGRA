
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')


### Load data
library(dplyr)

# RI data
load('9.2-Liver_RI_N=87.RData')

# BEQ data
load('9.1-Liver_XM_BEQ_CB_N=35.RData')

BEQ <- Liver_XM_BEQ_CB %>%
  rename(DTXSID = dtxsid) %>%
  rename_with(~ paste0("BEQ_", .), .cols = -1) #35 compounds
   
# BER data
load('9.4-Liver_HTTK_BER_N=47.RData')
BER <- Liver_BER %>%
  select(DTXSID, BER) #47 compounds

#### Merge RI 和 BEQ, BER ####
df <- Liver_RI %>%
  left_join(BER, by = 'DTXSID') %>%
  left_join(BEQ, by = 'DTXSID') 

unique(df$RI)
unique(df$BER)

#### ExpTox calculation ####
library(toxpiR)

# View all functions
lsf.str("package:toxpiR")

## Create slice with transformation functions
# Slice 1: BEQ, # log10(x)-log10(min)
# Slice 2: BER, # -log10(x)+log10(max)
# Slice 3: RI, # linear

# Automatically extract all column names starting with 'BEQ_'
BEQ_columns <- grep("^BEQ_", names(df), value = TRUE)

f.slices <- TxpSliceList(
  Slice1 = TxpSlice(BEQ_columns),
  Slice2 = TxpSlice('BER'),
  Slice3 = TxpSlice('RI')
)

## Method#1: Not consider missing value
linear <- function(x) {
  # Example linear transformation
  return(x)
}

final.trans <- TxpTransFuncList(
  f1 = function(x) {
    if (all(is.na(x))) return(x)
    log10(x) - log10(min(x, na.rm = TRUE))
  },
  f2 = function(x) {
    if (all(is.na(x))) return(x)
    -log10(x) + log10(max(x, na.rm = TRUE))
  },
  f3 = linear
)


## Goal - Create ToxPi model.
# All Slices, weight = 1.
f.model <- TxpModel(txpSlices = f.slices,
                    txpWeights = c (1,1,1),
                    txpTransFuncs = final.trans)


## Calculate ToxPi scores
f.results <- txpCalculateScores(model = f.model,
                                input = df,
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
# List all grobs in the current graph
grid.ls()
# Find all grobs containing 'label' and change the font size
all_grobs <- grid.ls(print = FALSE)


## Modify the font size of all labels
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
df_1 <- df %>%
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
Liver_ExpTox_Slices <- df_2 %>%
  full_join(toxpi_scores, by = 'id')

library(writexl)
write_xlsx(Liver_ExpTox_Slices, '10-Liver_ExpTox_Slice_N=94.xlsx')

#### Upgrade Rank plot ####
library(ggplot2)

# 去掉缺失值
Liver_ExpTox_Slices <- Liver_ExpTox_Slices %>%
  filter(!is.na(toxpi_score))

RankPlot <- ggplot(Liver_ExpTox_Slices, 
                   aes(x = toxpi_score, y = reorder(PREFERRED_NAME, toxpi_score))) +
  geom_point() +
  scale_x_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.1)) +
  labs(x = "ExpTox Score", y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(
      hjust = 1,  # 右对齐
      size = 6.0    # 减小字体大小，可以根据需要调整这个值
    )
  )

show(RankPlot)

# Output plot
ggsave("10.2-RankPlot_ExpTox_Liver.pdf", plot = RankPlot, 
       width = 6, height = 8, units = "in", dpi = 600)



### Hierarchical Clustering ####
f.hc <- hclust(dist(txpSliceScores(f.results)))
plot(f.hc, hang = -1, labels = txpIDs(f.results), xlab="", sub="")



