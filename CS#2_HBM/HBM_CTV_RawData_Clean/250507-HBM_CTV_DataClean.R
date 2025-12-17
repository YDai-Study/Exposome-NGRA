

setwd('/Users/ydai/Desktop/Exposome_NGRA/2505-HBM_Analysis/HBM_CTV_RawData_Clean/HBM_CTV_Pre')

# Load CTV raw data, whinch include previsous more compounds, also include new raw data
CTV_files <- list.files(pattern = "^ToxValue_Results_.*\\.csv$")

# 使用 lapply 读取所有CSV文件，返回一个包含数据框的列表
CTV_list <- lapply(CTV_files, read.csv)

# 如果需要将所有数据合并为一个数据框（假设列名相同）
CTV_merge <- do.call(rbind, CTV_list)

write.csv(CTV_merge, file = 'CTV_RawData_Pre_Merge.csv')

## Mannually merge CTV prediction and experiment data

# Load merged data
library(readxl)
CTV_Pred <- read_excel('/Users/ydai/Desktop/Exposome_NGRA/2505-HBM_Analysis/HBM_CTV_RawData/CTV_RawData_Exp_Pre_Merge.xlsx',
                  sheet = 'Prediction')

CTV_Exp <- read_excel('/Users/ydai/Desktop/Exposome_NGRA/2505-HBM_Analysis/HBM_CTV_RawData/CTV_RawData_Exp_Pre_Merge.xlsx',
                       sheet = 'Experiment')
CTV_Che <- read_excel('/Users/ydai/Desktop/Exposome_NGRA/2505-HBM_Analysis/HBM_CTV_RawData/CTV_RawData_Exp_Pre_Merge.xlsx',
                       sheet = 3)

### Clean CTV merged data
library(dplyr)
CTV_Che_1 <- CTV_Che %>%
  select(Chemical.name, DTXSID)

# Input DTXSID
CTV_Pred <- CTV_Pred %>%
  select(-...1)
CTV_Pred_1 <- left_join(CTV_Pred, CTV_Che_1, by = 'Chemical.name')

CTV_Exp <- CTV_Exp %>%
  dplyr::rename(Chemical.name = `Chemical name`)
CTV_Exp_1 <- left_join(CTV_Exp, CTV_Che_1, by = 'Chemical.name')

### Clean
CTV_Pred_2 <- CTV_Pred_1 %>%
  select(Chemical.name, DTXSID, `Model Name`, Unit, Prediction) %>%
  rename(Model_Name = `Model Name`) %>%
  filter(Unit %in% c('mg/(kg·day)', "mg/m3",
                     'risk per mg/(kg·day)', 'risk per µg/m3'))
unique(CTV_Pred_2$DTXSID) #65 compounds

CTV_Exp_2 <- CTV_Exp_1 %>%
  select(Chemical.name, DTXSID, Endpoint, `Toxicity value`, Unit) %>%
  rename(Toxicity_Value = `Toxicity value`) 
unique(CTV_Exp_2$DTXSID) #6 compounds

### Filter
HBM_76 <- read_excel('/Users/ydai/Desktop/Exposome_NGRA/2505-HBM_Analysis/250505-Overlap_EU_US_CA_N=76.xlsx',
                     sheet = 1) %>%
  select(DTXSID, PREFERRED_NAME)

CTV_Pred_3 <- left_join(HBM_76, CTV_Pred_2, by = 'DTXSID') %>%
  filter(!is.na(Prediction))
unique(CTV_Pred_3$DTXSID) #63 compounds

CTV_Exp_3 <- left_join(HBM_76, CTV_Exp_2, by = 'DTXSID') %>%
  filter(!is.na(Toxicity_Value))
unique(CTV_Exp_3$DTXSID) #5 compounds

### Save
setwd('/Users/ydai/Desktop/Exposome_NGRA/2505-HBM_Analysis/HBM_CTV_RawData_Clean')
CTV_Prediction <- CTV_Pred_3 
CTV_Experiment <- CTV_Exp_3
save(CTV_Prediction, CTV_Experiment, file = 'HBM_CTV_Pre_Exp.RData')