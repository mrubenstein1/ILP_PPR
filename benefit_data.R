
##### Benefit (Duck Density) Data ####
    # Import and clean benefit data
    # add to the established data panel


#### I. Import data ####
benefit <- read.csv("input_data/benefit.csv")

# restructure around gcm/rcp/time step structure
benefit_new <- benefit %>%
  mutate(
    year = case_when(
      scenario == "Obs" ~ 2020L,
      str_detect(scenario, "3059") ~ 2050L,  # mid‑century window
      str_detect(scenario, "7099") ~ 2100L,  # end‑century window
      TRUE ~ NA_integer_
    ),
    rcp = case_when(
      str_detect(scenario, "4.5") ~ "45",
      str_detect(scenario, "8.5") ~ "85",
      TRUE ~ NA_character_
    ),
    gcm = case_when(
      str_detect(str_to_lower(scenario), "dry") ~ "dry",
      str_detect(str_to_lower(scenario), "wet") ~ "wet",
      TRUE ~ NA_character_
    ),
    wmd_id = WMD,
    ben    = total_sum
  ) %>%
  select(wmd_id, year, rcp, gcm, ben)


#### II. Interpolate Decadal Data ####

library(dplyr)
library(tidyr)
library(purrr)

# 1. Start from benefit_new (with wmd_id, year, rcp, gcm, ben)
#    Create the decadal grid for each WMD × scenario
benefit_decades <- benefit_new %>%
  # Ensure year is numeric
  mutate(year = as.integer(year)) %>%
  # Create all decadal years 2020, 2030, ..., 2100 per group
  group_by(wmd_id, rcp, gcm) %>%
  complete(year = seq(2020L, 2100L, by = 10L)) %>%
  ungroup()

# 2. Linear interpolation of ben within each WMD × scenario group
benefit_interp <- benefit_decades %>%
  group_by(wmd_id, rcp, gcm) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    ben = {
      # x = observed years where ben is not NA
      x_obs  <- year[!is.na(ben)]
      y_obs  <- ben[!is.na(ben)]
      # if we have at least two points, interpolate; otherwise keep as is
      if (length(x_obs) >= 2) {
        approx(x = x_obs, y = y_obs, xout = year, method = "linear", rule = 2)$y
      } else {
        # no interpolation possible; return ben unchanged
        ben
      }
    }
  ) %>%
  ungroup()

# benefit_interp is your revised benefit_new with decadal, linearly interpolated ben