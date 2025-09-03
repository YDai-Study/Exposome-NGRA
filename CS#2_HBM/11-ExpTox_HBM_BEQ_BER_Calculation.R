
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

library(readxl)
library(dplyr)

#### API-CompTox-ctxR package, extract Target assay data ####

library(ctxR)
### API key will be stored as the variable my_key
my_key <- '5e944183-665f-4e1b-9ff0-c7e27b8b09e4'

### Retrieve all assays
assays <- get_all_assays(API_key = my_key) #1499 assays

library(writexl)
write_xlsx(assays, 'CompTox_assays.xlsx')

# Load Target Assays
TargetAssay <- read_excel('TargetedAssays.xlsx',
                          sheet = "HBM")

# Extract AEID for Target Assays and convert to Values
TA_AEID <- TargetAssay %>%
  select(aeid) %>%
  distinct() %>%
  pull(aeid) #12 AEIDs

# Extract DTXSID from HBM 
HBM <- read_excel('1-Overlap_EU_US_CA_N=76.xlsx',
                  sheet = 1)
DTXSID <- HBM %>%
  select(DTXSID) %>%
  distinct() %>%
  pull(DTXSID) #76 DTXSID

### Retrieve bioactivity data from NR_AEID batch
TA_bio_HBM <- get_bioactivity_details_batch(AEID = c(785,786,761,764,762,1816,
                                                     802,801,1198,1128,804,805), #12 Target Assays
                                            API_key = my_key,
                                            DTXSID = DTXSID #HBM compounds
)
save(TA_bio_HBM, file = "11.1-HBM_TargetAssay_Bioactivity.RData")

### Merge multiple lists into a database
library(dplyr)
library(purrr)

# Function to check if a dataframe is valid (non-empty)
is_valid_df <- function(df) {
  if (is.data.frame(df) && nrow(df) > 0 && ncol(df) > 0) {
    return(TRUE)
  }
  return(FALSE)
}

# Filter out empty dataframes from the list
load("11.1-HBM_TargetAssay_Bioactivity.RData")
TA_bio_HBM_filtered <- map(TA_bio_HBM, ~ {
  if (is_valid_df(.x)) .x else NULL
})

# Remove NULL elements from the list
TA_bio_HBM_filtered <- compact(TA_bio_HBM_filtered)

TA_bio_HBM_merge <- TA_bio_HBM_filtered %>%
  map_df(~ .x) #33480 records

unique(TA_bio_HBM_merge$dtxsid) #48 compounds

unique(TA_bio_HBM_merge$aeid)

# Save data
save(TA_bio_HBM_merge, file = '11.2-HBM_TA_Bioactivity_N=48.RData')

unique(TA_bio_HBM_merge$dtxsid) #48 compounds

##### Calculate BEQ ####
## Select target aeid 
load("11.2-HBM_TA_Bioactivity_N=48.RData")
HBM_TA_bio <- TA_bio_HBM_merge %>%
  semi_join(TargetAssay, by = 'aeid') #619 records

unique(HBM_TA_bio$aeid) #12 aeid
unique(HBM_TA_bio$dtxsid) #36 compounds

HBM_TA_bio_1 <- HBM_TA_bio %>%
  select(dtxsid, aeid, ac50)

## For records with the same dtxsid and aeid, select the record with the smallest ac50.
HBM_TA_bio_2 <- HBM_TA_bio_1 %>%
  group_by(dtxsid, aeid) %>%
  filter(ac50 ==  min(ac50)) %>%
  sample_n(1) %>%
  ungroup()  

## Calculate the minimum_ac50 for each assay as the AC50_ref
min_ac50 <- HBM_TA_bio %>%
  group_by(aeid) %>%
  summarize(ac50_ref = min(ac50, na.rm = TRUE))

df <- left_join(HBM_TA_bio_2, min_ac50, by = 'aeid')

#### Using CB to calculate BEQ ####
## Load HBM_CB data
load('7-HBM_Pred_Expo_CytoAC50_N=76.RData')
unique(HBM_Pred_Expo_CytoAC50$DTXSID)

HBM_CB <- HBM_Pred_Expo_CytoAC50 %>%
  select(DTXSID, CB_P50_uM) %>%
  rename(dtxsid = DTXSID,
         CB =  CB_P50_uM) %>%
  distinct()

# Input HBM_CB data
df.1 <- left_join(df, HBM_CB, by = 'dtxsid')
unique(df.1$dtxsid) #36 compounds

## Calculate CB/AC50 for each compound
CB_ac50 <- df.1 %>%
  group_by(dtxsid, aeid) %>%
  reframe(CB_ac50 = CB / ac50)

# Merge CB/AC50 and min_AC50
df.2 <- left_join(CB_ac50, min_ac50, by = 'aeid')

## Calculate BEQi = CBi/AC50i x min_AC50 for each compound
df.3 <- df.2 %>% 
  mutate(BEQ = CB_ac50 * ac50_ref)

unique(df.3$dtxsid) #36 compounds
unique(df.3$aeid) #12 assays

## Calculate BEQi% in each assay (aeid)
# Calculate sum_BEQ of all compounds in each assay
df.4 <- df.3 %>%
  group_by(aeid) %>%
  mutate(BEQ_sum = sum(BEQ, na.rm = TRUE)) %>%
  ungroup()

# Calculate BEQ% = (BEQi/BEQ_sum) x 100%
df.5 <- df.4 %>%
  mutate(BEQ_percent = (BEQ/BEQ_sum) * 100)

# Each compound conresponding to multiple target assays
library(tidyr)
df.6 <- df.5 %>%
  select(dtxsid, aeid, BEQ_percent) %>%
  pivot_wider(names_from = aeid,
              values_from = BEQ_percent)

HBM_TA_BEQ_CB <- df.6
save(HBM_TA_BEQ_CB, file = '11.3-HBM_TA_BEQ_CB_N=36.RData')

#### Using Urine_Conc to calculate BEQ ####

## Load HBM_Conc data
load('7-HBM_Pred_Expo_CytoAC50_N=76.RData')

HBM <- HBM_Pred_Expo_CytoAC50 %>%
  select(InChIKey, Sample_Detail, Chemical_Name, DTXSID, Chemical_Group,
         P5_Conc_P50_uM, P50_Conc_P50_uM, P95_Conc_P50_uM) %>%
  filter(Sample_Detail == 'Urine') %>%
  filter(Chemical_Group != 'metals and trace elements') %>%
  distinct() %>%
  filter(!is.na(P50_Conc_P50_uM))

HBM_Urine <- HBM %>%
  select(DTXSID, P50_Conc_P50_uM) %>%
  rename(dtxsid = DTXSID) #46 compounds
unique(HBM_Urine$dtxsid)

# Input HBM_Urine data
df.1 <- left_join(df, HBM_Urine, by = 'dtxsid')
unique(df.1$dtxsid) #36 compounds

## Calculate Conc_uM_P50/AC50 for each compound
Conc_ac50 <- df.1 %>%
  group_by(dtxsid, aeid) %>%
  reframe(Conc_ac50 = P50_Conc_P50_uM / ac50)

# Merge Conc_uM_P50/AC50 and min_AC50
df.2 <- left_join(Conc_ac50, min_ac50, by = 'aeid')

## Calculate BEQi = Conci/AC50i x min_AC50 for each compound
df.3 <- df.2 %>% 
  mutate(BEQ = Conc_ac50 * ac50_ref)

unique(df.3$dtxsid) #36 compounds
unique(df.3$aeid) #12 assays

## Calculate BEQi% in each assay (aeid)
# Calculate sum_BEQ of all compounds in each assay
df.4 <- df.3 %>%
  group_by(aeid) %>%
  mutate(BEQ_sum = sum(BEQ, na.rm = TRUE)) %>%
  ungroup()

# Calculate BEQ% = (BEQi/BEQ_sum) x 100%
df.5 <- df.4 %>%
  mutate(BEQ_percent = (BEQ/BEQ_sum) * 100) %>%
  filter(!is.na(BEQ_percent))

unique(df.5$dtxsid) #27 compounds

# Each compound conresponding to multiple target assays
library(tidyr)
df.6 <- df.5 %>%
  select(dtxsid, aeid, BEQ_percent) %>%
  pivot_wider(names_from = aeid,
              values_from = BEQ_percent)

HBM_TA_BEQ_Urine <- df.6
save(HBM_TA_BEQ_Urine, file = '11.4-HBM_TA_BEQ_Urine_N=27.RData')


##### Calculate BER #####

## Load HBM Overlapping compounds N=76 ##
library(readxl)
Che_list <- read_excel('1-Overlap_EU_US_CA_N=76.xlsx',
                       sheet = 1)

### HTTK model to calculate Css ###
library(httk)
library(readxl)
library(ggplot2)

### Step 1: Loading In Vitro Screening Data ####

load('invitrodb_3_5_mc5.Rdata')

# Built new subset of chemicals with human HTTK data
toxcast.httk <- subset(mc5, dsstox_substance_id %in% get_cheminfo(
  info='DTXSID', suppress.messages = TRUE))

# Extract 'DTXSID' from our chemcial list
my.chems <- Che_list$DTXSID #76

# Extract cooresponding chemicals from Toxcast
example.toxcast <- as.data.frame(mc5[mc5$dsstox_substance_id %in% my.chems,]) #51074 records

# pick 9 variables
example.toxcast <- example.toxcast[, c("chnm",
                                       "dsstox_substance_id",
                                       "spid",
                                       "hitc",
                                       "modl",
                                       "aeid",
                                       "modl_ga",
                                       "modl_ac10",
                                       "modl_acc")]

# reduce precision to decrease space:
cols <- c("modl_ga", "modl_ac10", "modl_acc")
for (this.col in cols)
  example.toxcast[, this.col] = signif(example.toxcast[, this.col], 3)

#Each chemical might have multiple assays for which there is a hit. 
#Each chemical may have been run multiple times (experimental replicates), 
#so that there may be multiple samples from the same chemical for the same assay. 
#Different samples are indicated by the column “spid” (sample id).

knitr::kable(head(example.toxcast), caption = "ToxCast In Vitro Bioactivity Data",
             floating.environment="sidewaystable")

# Find the in vitro concentrations for compounds
toxcast.table <- NULL 
for (this.id in unique(example.toxcast$dsstox_substance_id))
{
  this.subset <- subset(example.toxcast, dsstox_substance_id == this.id)
  these.hits <- subset(this.subset, hitc==1)
  if (dim(these.hits)[1]>0){
    this.row <- data.frame(Compound=as.character(this.subset[1,"chnm"]),
                           DTXSID=this.id,
                           Total.Assays = dim(this.subset)[1],
                           Unique.Assays = length(unique(this.subset$aeid)),
                           Total.Hits = dim(these.hits)[1],
                           Unique.Hits = length(unique(these.hits$aeid)),
                           Low.AC50 = signif(min(these.hits$modl_ga),3),
                           Low.AC10 = signif(min(these.hits$modl_ac10),3),
                           Low.ACC = signif(min(these.hits$modl_acc),3),
                           Q10.AC50 = signif(quantile(these.hits$modl_ga,probs=0.1),3),
                           Q10.AC10 = signif(quantile(these.hits$modl_ac10,probs=0.1),3),
                           Q10.ACC = signif(quantile(these.hits$modl_acc,probs=0.1),3),
                           Med.AC50 = signif(median(these.hits$modl_ga),3),
                           Med.AC10 = signif(median(these.hits$modl_ac10),3),
                           Med.ACC = signif(median(these.hits$modl_acc),3),
                           Q90.AC50 = signif(quantile(these.hits$modl_ga,probs=0.9),3),
                           Q90.AC10 = signif(quantile(these.hits$modl_ac10,probs=0.9),3),
                           Q90.ACC = signif(quantile(these.hits$modl_acc,probs=0.9),3)
    )
    toxcast.table <- rbind(toxcast.table, this.row)
  }
}
rownames(toxcast.table) <- seq(1,dim(toxcast.table)[1]) #39
knitr::kable(head(toxcast.table[,1:6]), caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")


# Save HBM_HTTK_toxcast.table
library(writexl)
write_xlsx(toxcast.table, '11-HBM_HTTK_toxcast_table.xlsx')

### Step 2: HTTK deduce Css_Plasma with in vitro ####

#Steady State Plasma Concentration Css 
#calc_mc_css(): predict Css resulting from an ongoing 1 mg/kg/day exposure.

# Assume toxcast.table already has columns "Css_95th", "Css_50th", "Css_5th"
for (this.id in unique(toxcast.table$DTXSID)) {
  # Check if a chemical substance is in the HTTK database
  if (this.id %in% get_cheminfo(info = "dtxsid", suppress.messages = TRUE)) {
    # Set the random number seed to ensure reproducibility of Monte Carlo simulations
    set.seed(12345)
    
    # Calculating CSS at different percentiles
    css_results <- calc_mc_css(dtxsid = this.id,
                               output.units = "uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # httk 默认采样数
                               suppress.messages = TRUE)
    
    # Store the results in a table
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_5th"]  <- css_results[3]
    
    # Mark CSS Type
    toxcast.table[toxcast.table$DTXSID == this.id, "Css.Type"] <- "in vitro"
  }
}

#Let’s look at just the relevant columns from our table. 
#Any Css values of “NA” (not a number) indicate that we don’t currently have in vitro TK data for those chemicals:

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th', "Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")


### Step 3: HTTK deduce Css_Plasma with QSPR Predictions ####

#For chemicals where no in vitro toxicokinetic data are available, 
#we can sometimes load predictions from in silico models. 

load_sipes2017()

library(httk)

# Make sure the Css column exists, if not add it
if (!"Css_95th" %in% colnames(toxcast.table)) {
  toxcast.table$Css_95th <- NA
}
if (!"Css_50th" %in% colnames(toxcast.table)) {
  toxcast.table$Css_50th <- NA
}
if (!"Css_5th" %in% colnames(toxcast.table)) {
  toxcast.table$Css_5th <- NA
}
if (!"Css.Type" %in% colnames(toxcast.table)) {
  toxcast.table$Css.Type <- NA
}

# Iterate over each chemical ID
for (this.id in unique(toxcast.table$DTXSID)) {
  
  # Check if the chemical is in the HTTK database and Css is empty
  if (this.id %in% get_cheminfo(info="dtxsid", suppress.messages=TRUE) &
      is.na(toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"])) {
    
    # 
    set.seed(12345)
    
    # Calculate 95th, 50th and 5th percentile CSS
    css_results <- calc_mc_css(dtxsid=this.id,
                               output.units="uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # 默认样本数
                               suppress.messages=TRUE)
    
    # Store calculation results in different columns
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_5th"] <- css_results[3]
    
    # Mark CSS type as "in silico"
    toxcast.table[toxcast.table$DTXSID==this.id, "Css.Type"] <- "in silico"
  }
}

#Take another look at our table – many of the gaps are now filled in.

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th',"Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")

sum(!is.na(toxcast.table$Css_50th))

## Output
HBM_Css_Plsama = toxcast.table

## Step 4: Reverse Dosimetry In Vitro-In Vivo Extropolation ####

#Don’t forget the ToxCast concentrations are given on the log-10 scale:
HBM_Css_Plsama$EquivDose <- signif(10^HBM_Css_Plsama$Q10.AC50 / HBM_Css_Plsama$Css_50th,
                                   3)

knitr::kable(HBM_Css_Plsama[1:10,c("Compound","Q10.AC50","Css_50th","EquivDose")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")          

## Step 5: Compare with Estimate Exposure Rates ####

# Load SEEM3 database
load('/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/External_Prediction.RData')

library(dplyr)
SEEM <- External_Prediction %>%
  select(DTXSID, InChIKey, P5, P50, P95)

# Pick interested chemicals
example.seem <- as.data.frame(subset(SEEM,
                                     DTXSID %in% HBM_Css_Plsama$DTXSID))

# Merge database
HBM_Css_Plsama <- left_join(HBM_Css_Plsama, example.seem, by = 'DTXSID')

#Calculate a bioactivity:exposure ratio (BER)
HBM_Css_Plsama$BER <- signif(HBM_Css_Plsama$EquivDose/HBM_Css_Plsama$P95,3)

HBM_BER <- HBM_Css_Plsama %>%
  select(InChIKey, DTXSID, BER) %>%
  filter(!is.na(BER))

unique(HBM_BER$InChIKey)
unique(HBM_BER$DTXSID)

save(HBM_BER, file = '11.5-HBM_HTTK_BER_N=20.RData')


