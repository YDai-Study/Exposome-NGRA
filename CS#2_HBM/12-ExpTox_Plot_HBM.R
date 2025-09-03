
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

### BEQ_CB TargetAssay ####
load('11.3-HBM_TA_BEQ_CB_N=36.RData')
BEQ <- HBM_TA_BEQ_CB %>%
  rename(DTXSID = dtxsid) %>%
  rename_with(~ paste0("BEQ_", .), .cols = -1) %>%
  # remove compounds with many NA values
  filter(!(DTXSID %in% c('DTXSID7020508', 'DTXSID9020376'))) #34 compounds

### BER ####
load('11.5-HBM_HTTK_BER_N=20.RData')
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




