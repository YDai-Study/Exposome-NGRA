
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#2_HBM')

## Load HBM Overlapping compounds N=76 ##
library(readxl)
Che_list <- read_excel('1-Overlap_EU_US_CA_N=76.xlsx',
                     sheet = 1)

### HTTK model to calculate Css ###
library(httk)
library(readxl)
library(ggplot2)

### Step 1: Loading In Vitro Screening Data ####

load('/Users/ydai/Library/CloudStorage/OneDrive-Personal/2-2023-UMCU-Postdoc/1-UMCU-Postdoc Project/4-Prediction Tools/HTTK_R_Code/invitrodb_3_5_mc5.Rdata')

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
write_xlsx(toxcast.table, '5.1-HBM_HTTK_toxcast_table.xlsx')

### Step 2: HTTK deduce Css_Plasma with in vitro ####

#Steady State Plasma Concentration Css 
#calc_mc_css(): predict Css resulting from an ongoing 1 mg/kg/day exposure.

# 假设 toxcast.table 中已经有列 "Css_95th", "Css_50th", "Css_5th"
for (this.id in unique(toxcast.table$DTXSID)) {
  # 检查化学物质是否在 HTTK 数据库中
  if (this.id %in% get_cheminfo(info = "dtxsid", suppress.messages = TRUE)) {
    # 设置随机数种子，以确保蒙特卡洛模拟的可重复性
    set.seed(12345)
    
    # 计算不同百分位数的 Css
    css_results <- calc_mc_css(dtxsid = this.id,
                               output.units = "uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # httk 默认采样数
                               suppress.messages = TRUE)
    
    # 将结果存入表格
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_5th"]  <- css_results[3]
    
    # 标记 Css 类型
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

# 遍历每个化学物质 ID
for (this.id in unique(toxcast.table$DTXSID)) {
  
  # 检查化学物质是否在 HTTK 数据库中，并且 Css 为空
  if (this.id %in% get_cheminfo(info="dtxsid", suppress.messages=TRUE) &
      is.na(toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"])) {
    
    # 设置随机数种子以保证结果可重复
    set.seed(12345)
    
    # 计算 95th, 50th 和 5th 百分位数的 Css
    css_results <- calc_mc_css(dtxsid=this.id,
                               output.units="uM",
                               concentration = 'plasma',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # 默认样本数
                               suppress.messages=TRUE)
    
    # 将计算结果存储到不同的列中
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_5th"] <- css_results[3]
    
    # 标记 Css 类型为 "in silico"
    toxcast.table[toxcast.table$DTXSID==this.id, "Css.Type"] <- "in silico"
  }
}

#Take another look at our table – many of the gaps are now filled in.

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th',"Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")

sum(!is.na(toxcast.table$Css_50th))

## Output
HBM_Css_Plasma = toxcast.table
HBM_Css_Plasma <- HBM_Css_Plasma %>%
  select(DTXSID, Css_95th, Css_50th, Css_5th, Css.Type) %>%
  dplyr::rename(Css_plasma_95th = Css_95th,
                Css_plasma_50th = Css_50th,
                Css_plasma_5th = Css_5th,
                Css.Type_Plasma = Css.Type)

### Deduce Blood Css ####

### Step 2.1: HTTK deduce Css_Blood with in vitro ####

#Steady State Plasma Concentration Css 
#calc_mc_css(): predict Css resulting from an ongoing 1 mg/kg/day exposure.

# 假设 toxcast.table 中已经有列 "Css_95th", "Css_50th", "Css_5th"
for (this.id in unique(toxcast.table$DTXSID)) {
  # 检查化学物质是否在 HTTK 数据库中
  if (this.id %in% get_cheminfo(info = "dtxsid", suppress.messages = TRUE)) {
    # 设置随机数种子，以确保蒙特卡洛模拟的可重复性
    set.seed(12345)
    
    # 计算不同百分位数的 Css
    css_results <- calc_mc_css(dtxsid = this.id,
                               output.units = "uM",
                               concentration = 'blood',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # httk 默认采样数
                               suppress.messages = TRUE)
    
    # 将结果存入表格
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID == this.id, "Css_5th"]  <- css_results[3]
    
    # 标记 Css 类型
    toxcast.table[toxcast.table$DTXSID == this.id, "Css.Type"] <- "in vitro"
  }
}

#Let’s look at just the relevant columns from our table. 
#Any Css values of “NA” (not a number) indicate that we don’t currently have in vitro TK data for those chemicals:

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th', "Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")


### Step 3.1: HTTK deduce Css_Blood with QSPR Predictions ####

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

# 遍历每个化学物质 ID
for (this.id in unique(toxcast.table$DTXSID)) {
  
  # 检查化学物质是否在 HTTK 数据库中，并且 Css 为空
  if (this.id %in% get_cheminfo(info="dtxsid", suppress.messages=TRUE) &
      is.na(toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"])) {
    
    # 设置随机数种子以保证结果可重复
    set.seed(12345)
    
    # 计算 95th, 50th 和 5th 百分位数的 Css
    css_results <- calc_mc_css(dtxsid=this.id,
                               output.units="uM",
                               concentration = 'blood',
                               which.quantile = c(0.95, 0.5, 0.05),
                               samples = 1000,  # 默认样本数
                               suppress.messages=TRUE)
    
    # 将计算结果存储到不同的列中
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_95th"] <- css_results[1]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_50th"] <- css_results[2]
    toxcast.table[toxcast.table$DTXSID==this.id, "Css_5th"] <- css_results[3]
    
    # 标记 Css 类型为 "in silico"
    toxcast.table[toxcast.table$DTXSID==this.id, "Css.Type"] <- "in silico"
  }
}

#Take another look at our table – many of the gaps are now filled in.

knitr::kable(toxcast.table[1:10,c("Compound","Q10.AC50","Css_95th",'Css_50th', 'Css_5th',"Css.Type")], 
             caption = "Summarized ToxCast Data",
             floating.environment="sidewaystable")

sum(!is.na(toxcast.table$Css_50th))

## Output
HBM_Css_Blood = toxcast.table
HBM_Css_Blood <- HBM_Css_Blood %>%
  select(DTXSID, Css_95th, Css_50th, Css_5th, Css.Type) %>%
  dplyr::rename(Css_blood_95th = Css_95th,
                Css_blood_50th = Css_50th,
                Css_blood_5th = Css_5th,
                Css.Type_Blood = Css.Type)

##### Merge plasma and blood Css ####
HBM_Css <- left_join(HBM_Css_Blood, HBM_Css_Plasma, by = 'DTXSID')

### Extract predicted internal blood from HExpPredict ####
load("/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/Internal_Prediction.RData")
Predict_Blood <- Internal_Prediction %>%
  select(DTXSID, P5, P50, P95) %>%
  dplyr::rename(CB_P5_uM = P5,
                CB_P50_uM = P50,
                CB_P95_uM = P95) 

HBM_Predict_Blood <- left_join(Che_list, Predict_Blood, by = 'DTXSID')

### Extract predicted external exposure from SEEM3 ####
load('/Users/ydai/Documents/GitHub/Exposome-NGRA/Exposome_Database/Clean_Data/External_Prediction.RData')
Predict_External <- External_Prediction %>%
  select(DTXSID, P5, P50, P95, Concentration_Unit) %>%
  dplyr::rename(Ex_P5 = P5,
                Ex_P50 = P50,
                Ex_P95 = P95,
                Ex_Unit = Concentration_Unit)

HBM_Predict <- left_join(HBM_Predict_Blood, Predict_External, 
                                  by = 'DTXSID')

##### Merge ####
Pred_In_Ex <- left_join(HBM_Predict, HBM_Css_Blood, by = 'DTXSID')
Pred_In_Ex_1 <- left_join(Pred_In_Ex, HBM_Css_Plsama, by = 'DTXSID')

library(writexl)
write_xlsx(Pred_In_Ex_1, '5.2-Overlap_HBM_Prediction.xlsx')