
setwd('/Users/ydai/Documents/GitHub/Exposome-NGRA/CS#1_Liver')
library(readxl)
library(dplyr)

Tox <- read_excel('8-Liver_Tox_Slices_N=94.xlsx',
                  sheet = 1)
Tox_1 <- Tox %>%
  select(PREFERRED_NAME, toxpi_score) %>%
  mutate(rank_tox = rank(-toxpi_score, ties.method = "average")) %>%
  rename(Tox_score = toxpi_score)

ExpTox <- read_excel('10-Liver_ExpTox_Slice_N=94.xlsx',
                     sheet = 1)
ExpTox_1 <- ExpTox %>%
  select(PREFERRED_NAME, toxpi_score) %>%
  mutate(rank_exp = rank(-toxpi_score, ties.method = 'average')) %>%
  rename(ExpTox_score = toxpi_score)

Tox_Exp <- left_join(Tox_1, ExpTox_1, by = 'PREFERRED_NAME')

Tox_Exp <- Tox_Exp %>%
  mutate(change = rank_tox - rank_exp)

# Plotting ranking changes

library(ggplot2)
library(ggrepel)

ggplot(Tox_Exp, aes(x = rank_tox, y = rank_exp, label = PREFERRED_NAME)) +
  geom_point(aes(color = change), size = 1.5) +
  geom_text_repel(size = 3.5, max.overlaps = Inf) + #  geom_text_repel
  scale_color_gradient2(low = "blue", high = "red", mid = "purple", midpoint = 0) +
  labs(title = "Rank Change Plot", x = "Tox_rank", y = "ExpTox_rank") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),      # 
    axis.title = element_text(size = 12),                    # 
    axis.text = element_text(size = 12),                     # 
    legend.title = element_text(size = 11),                  # 
    legend.text = element_text(size = 11)                    # 
  )

