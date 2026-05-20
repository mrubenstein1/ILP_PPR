library(dplyr)


#impose 10% transition rate uniform
eau_panel_alloc <- eau_panel_alloc %>%
  mutate(trans_prob = 0.10)