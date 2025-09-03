
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

#### API-CompTox-ctxR package, extract Liver_Xenobiotic metabolism data ####

library(ctxR)
### API key will be stored as the variable my_key
my_key <- '5e944183-665f-4e1b-9ff0-c7e27b8b09e4'

### Retrieve all assays
assays <- get_all_assays(API_key = my_key) #1499 assays

## Select human, liver, nuclear receptor
library(dplyr)
assays_1 <- assays %>%
  filter(organism == 'human' &
           tissue == 'liver' &
           intendedTargetFamilySub == 'xenobiotic metabolism') #34 assays

XM_aeid <- assays_1 %>%
  select(aeid)

library(writexl)
write_xlsx(assays_1, '9.1-Liver_XenoMetabolism_assays.xlsx')

# Extract DTXSID from Liver 
Liver <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                    sheet = 1)
DTXSID <- Liver %>%
  select(DTXSID) %>%
  distinct() %>%
  pull(DTXSID) #94 DTXSID

### Retrieve bioactivity data from NR_AEID batch
XM_bio_Liver <- get_bioactivity_details_batch(
  API_key = my_key,
  DTXSID = DTXSID #Liver compounds
)


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
XM_bio_Liver_filtered <- map(XM_bio_Liver, ~ {
  if (is_valid_df(.x)) .x else NULL
})

# Remove NULL elements from the list
XM_bio_Liver_filtered <- compact(XM_bio_Liver_filtered)

XM_bio_Liver_merge <- XM_bio_Liver_filtered %>%
  map_df(~ .x) #49215 records

# XM_aeid_bio_Liver
XM_bio_Liver_merge_1 <- XM_bio_Liver_merge %>%
  select(dtxsid, aeid, ac50)

XM_aeid_bio_Liver <- left_join(XM_aeid, XM_bio_Liver_merge_1, by = 'aeid')


unique(XM_aeid_bio_Liver$dtxsid) #35 compounds

unique(XM_aeid_bio_Liver$aeid) #34 XM assays

# Save data
library(writexl)
write_xlsx(XM_aeid_bio_Liver, '9.2-XM_aeid_bio_Liver_N=35.xlsx') 

#### BEQ calcualtion ####
library(readxl)
Liver_XM_bio <- read_excel('9.2-XM_aeid_bio_Liver_N=35.xlsx',
                           sheet = 1) #548 records

## For records with the same dtxsid and aeid, select the record with the smallest ac50
Liver_XM_bio_1 <- Liver_XM_bio %>%
  group_by(dtxsid, aeid) %>%
  filter(ac50 ==  min(ac50)) %>%
  sample_n(1) %>%
  ungroup()  # 取消分组

## Calculate the minimum_ac50 for each assay as the AC50_ref
min_ac50 <- Liver_XM_bio %>%
  group_by(aeid) %>%
  summarize(ac50_ref = if (any(!is.na(ac50))) min(ac50, na.rm = TRUE) else NA_real_)


df <- left_join(Liver_XM_bio_1, min_ac50, by = 'aeid')

#### Using CB to calculate BEQ ####
## Load Liver_CB data
load('5-Liver_Pred_Expo_CytoAC50_N=94.RData')

Liver_CB <- Liver_Pred_Expo_CytoAC50 %>%
  select(DTXSID,  CB_P50_uM ) %>%
  distinct() %>%
  filter(!is.na(CB_P50_uM)) %>% 
  rename(dtxsid = DTXSID) #78 compounds

# Input HBM_CB data
df.1 <- left_join(df, Liver_CB, by = 'dtxsid')
unique(df.1$dtxsid) #35 compounds

## Calculate CB/AC50 for each compound
CB_ac50 <- df.1 %>%
  group_by(dtxsid, aeid) %>%
  reframe(CB_ac50 = CB_P50_uM / ac50)

# Merge CB/AC50 and min_AC50
df.2 <- left_join(CB_ac50, min_ac50, by = 'aeid')

## Calculate BEQi = CBi/AC50i x min_AC50 for each compound
df.3 <- df.2 %>% 
  mutate(BEQ = CB_ac50 * ac50_ref)

unique(df.3$dtxsid) #35 compounds
unique(df.3$aeid) #32 assays

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

df.7 <- df.6 %>%
  select(1:(which(names(df.6) == "1398") - 1)) #15 assays

Liver_XM_BEQ_CB <- df.7
save(Liver_XM_BEQ_CB, file = '9.1-Liver_XM_BEQ_CB_N=35.RData')


### RI calculation using SEEM3 ####
library(readxl)
Toxpi <- read_excel('8-Liver_Tox_Slices_N=94.xlsx',
                    sheet = 1)

load('5-Liver_Pred_Expo_CytoAC50_N=94.RData')

library(dplyr)
Toxpi <- Toxpi %>%
  select(DTXSID, PREFERRED_NAME, toxpi_score)

SEEM3 <- Liver_Pred_Expo_CytoAC50 %>%
  select(DTXSID, Ex_P50) %>%
  distinct() %>%
  filter(!is.na(Ex_P50)) %>%
  mutate(log_seem3 = log(Ex_P50)) # log-transformed SEEM3

df <- left_join(Toxpi, SEEM3, by = 'DTXSID')

## Normalized SEEM3

df$seem3_normalized <- (df$log_seem3 - min(df$log_seem3, na.rm = TRUE)) / 
  (max(df$log_seem3, na.rm = TRUE) - min(df$log_seem3, na.rm = TRUE))

## Calulate RI
df <- df %>%
  mutate(RI = seem3_normalized * toxpi_score) #87 compounds

# Save data
Liver_RI <- df
save(Liver_RI, file = '9.2-Liver_RI_N=87.RData')

##### BER Calculation using HTTK #####
library(readxl)

Che_list <- read_excel('1.1-Overlap_IARC_NTP_N=94.xlsx',
                       sheet = 1)

### HTTK model to calculate Css ###
library(httk)
library(readxl)
library(ggplot2)
library(dplyr)

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

### Step 2: HTTK deduce Css_Plasma with in vitro ####

#Steady State Plasma Concentration Css 
#calc_mc_css(): predict Css resulting from an ongoing 1 mg/kg/day exposure.

# Assume toxcast.table already has columns "Css_95th", "Css_50th", "Css_5th"
for (this.id in unique(toxcast.table$DTXSID)) {
  # Check if a chemical substance is in the HTTK database
  if (this.id %in% get_cheminfo(info = "dtxsid", suppress.messages = TRUE)) {
    # 
    set.seed(12345)
    
    # Calculating Css at different percentiles
    css_results <- calc_mc_css(dtxsid = this.id,
                               output.units = "uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # 
                               suppress.messages = TRUE)
    
    # Store the results in a table
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_5th"]  <- css_results[3]
    
    # Tag Css Type
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

# 确保 Css 列存在，如果没有则添加
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
    
    # Set the random number seed to ensure repeatable results
    set.seed(12345)
    
    # Calculate 95th, 50th and 5th percentile Css
    css_results <- calc_mc_css(dtxsid=this.id,
                               output.units="uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # Default number of samples
                               suppress.messages=TRUE)
    
    # Store calculation results in different columns
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_5th"] <- css_results[3]
    
    # Mark Css type as "in silico"
    toxcast.table[toxcast.table$DTXSID==this.id, "Css.Type"] <- "in silico"
  }
}

#Take another look at our table – many of the gaps are now filled in.

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th',"Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")

sum(!is.na(toxcast.table$Css_50th))

## Output
Liver_Css_Plasma = toxcast.table %>%
  filter(!is.na(Css_50th)) #47 compounds

save(Liver_Css_Plasma, file = '9.3-Liver_Css_Plasama_N=47.RData')

## Step 4: Reverse Dosimetry In Vitro-In Vivo Extropolation ####
load('9.3-Liver_Css_Plasama_N=47.RData')
#Don’t forget the ToxCast concentrations are given on the log-10 scale:
Liver_Css_Plasma$EquivDose <- signif(10^Liver_Css_Plasma$Q10.AC50 / Liver_Css_Plasma$Css_50th,
                                     3)

knitr::kable(Liver_Css_Plasma[1:10,c("Compound","Q10.AC50","Css_50th","EquivDose")], 
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
                                     DTXSID %in% Liver_Css_Plasma$DTXSID))

# Merge database
Liver_Css_Plasma <- left_join(Liver_Css_Plasma, example.seem, by = 'DTXSID')

#Calculate a bioactivity:exposure ratio (BER)
Liver_Css_Plasma$BER <- signif(Liver_Css_Plasma$EquivDose/Liver_Css_Plasma$P95,3)

Liver_BER <- Liver_Css_Plasma %>%
  select(InChIKey, DTXSID, BER) %>%
  filter(!is.na(BER))

unique(Liver_BER$InChIKey)
unique(Liver_BER$DTXSID) #47 compounds

save(Liver_BER, file = '9.4-Liver_HTTK_BER_N=47.RData')


