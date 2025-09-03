
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')

library(readxl)
library(dplyr)
df <- read_excel('5-Overlap_Liver_Conc_uM_Unit_P5.P50.P95.xlsx',
                         sheet = 1)
unique(df$InChIKey) #94 compounds

df_1 <- df %>%
  mutate(
    Detection = case_when(
      !is.na(Sample) ~ "Yes",
      TRUE ~ "No"
    ),
    Available_Conc_uM = case_when(
      !is.na(Sample) & !is.na(Conc_P50_uM) ~ "Yes",
      TRUE ~ "No"
    )
  )

df_2 <- df_1 %>%
  select(InChIKey, Sample, Detection, Available_Conc_uM) %>%
  mutate(
    Available_Conc_uM_1 = case_when(
      Detection == "Yes" ~ paste0(Available_Conc_uM, ",", Sample),
      TRUE ~ Available_Conc_uM
    )
  )
unique(df_2$InChIKey)

library(stringr)

df_2 <- df_1 %>%
  select(InChIKey, Sample, Detection, Available_Conc_uM) %>%
  group_by(InChIKey) %>%
  mutate(
    Sample_1 = if(any(Detection == "Yes")) {
      str_c(unique(Sample), collapse = ",")
    } else {
      NA_character_  
    }
  ) %>%
  ungroup()

df_3 <- df_2 %>%
  group_by(InChIKey) %>%
  mutate(
    Available_Conc_uM_1 = case_when(
      Available_Conc_uM == "Yes" ~ {
        
        samples_concat <- str_c(unique(Sample[Available_Conc_uM == "Yes"]), collapse = ",")
     
        str_c("Yes", samples_concat, sep = ",")
      },
      Available_Conc_uM == "No" ~ NA_character_,
      TRUE ~ NA_character_  
    )
  ) %>%
  ungroup()

df_4 <- df_3 %>%
  select(-Sample, -Available_Conc_uM) %>%
  distinct()

unique(df_4$InChIKey) #94 compounds


df_5 <- df_4 %>%
  group_by(InChIKey, Sample_1) %>%    
  filter(
  
    !(Detection == "Yes" & n() > 1 & is.na(Available_Conc_uM_1))
  ) %>%
  ungroup() %>%
  rename(Sample = Sample_1, Available_Conc_uM = Available_Conc_uM_1)

library(writexl)
write_xlsx(df_5, '5.1-Overlap_Liver_Sample.xlsx')

