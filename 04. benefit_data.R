
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


# 0. Keep track of which rows are observed vs modeled future
benefit_tagged <- benefit_new %>%
  mutate(
    is_obs = is.na(rcp) & is.na(gcm) & year == 2020L
  )

# 1. Split into observed baseline and future scenarios
benefit_obs <- benefit_tagged %>%
  filter(is_obs)

benefit_future <- benefit_tagged %>%
  filter(!is_obs)  # only rows with RCP/GCM defined

# 2. For future scenarios, create decadal grid and interpolate within each WMD × scenario

benefit_future_decades <- benefit_future %>%
  mutate(year = as.integer(year)) %>%
  group_by(wmd_id, rcp, gcm) %>%
  complete(year = seq(2030L, 2100L, by = 10L)) %>%   # full decadal grid
  ungroup()

benefit_future_interp <- benefit_future_decades %>%
  group_by(wmd_id, rcp, gcm) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    ben = {
      x_obs <- year[!is.na(ben)]
      y_obs <- ben[!is.na(ben)]
      if (length(x_obs) >= 2) {
        # Linear interpolation between the available anchors
        approx(x = x_obs, y = y_obs, xout = year,
               method = "linear", rule = 2)$y
      } else {
        # If only one point exists (e.g., only 2050), we leave others as NA
        ben
      }
    }
  ) %>%
  ungroup()

# 3. Recombine observed baseline and interpolated future scenarios

benefit_interp <- bind_rows(
  benefit_obs,          # only 2020 Obs values
  benefit_future_interp # decadal interpolated futures
) %>%
  arrange(wmd_id, rcp, gcm, year)



#####III. Match Benefit data into EAU_Panel ####


# 1. Select only the join keys + ben from benefit_interp
benefit_for_join <- benefit_interp %>%
  select(wmd_id, year, rcp, gcm, ben)

# 2. Join to eau_panel and write ben into abs_density
eau_panel_with_ben <- eau_panel %>%
  left_join(benefit_for_join,
            by = c("wmd_id", "year", "rcp", "gcm")) %>%
  mutate(
    abs_density = ben  # copy WMD-level ben into abs_density
  ) %>%
  select(-ben)  # optional: drop ben after copying

### check that join worked properly

# For each WMD × year × scenario, check if all EAUs share the same abs_density
abs_check <- eau_panel_with_ben %>%
  group_by(wmd_id, year, rcp, gcm) %>%
  summarise(
    n_eau        = n(),
    min_abs      = min(abs_density, na.rm = TRUE),
    max_abs      = max(abs_density, na.rm = TRUE),
    n_na_abs     = sum(is.na(abs_density)),
    .groups = "drop"
  ) %>%
  mutate(
    all_equal = (min_abs == max_abs) | n_eau <= 1,  # or only one EAU
    any_na    = n_na_abs > 0
  )

# Look at any combinations where abs_density is not constant across EAUs
abs_check %>%
  filter(!all_equal & !any_na)

####. Distribute benefit by proportional habitat #####

eau_panel_alloc <- eau_panel_with_ben %>%
  group_by(wmd_id, year, rcp, gcm) %>%
  mutate(
    # total suitable habitat in this WMD × year × scenario
    total_prop_suitable = sum(prop_suitable, na.rm = TRUE),
    # share of WMD habitat carried by this EAU
    habitat_share = ifelse(
      total_prop_suitable > 0,
      prop_suitable / total_prop_suitable,
      NA_real_
    ),
    # allocate WMD-level abs_density down to EAUs
    scaled_density = habitat_share * abs_density
  ) %>%
  ungroup()

### check that this distribution worked propoerly 
alloc_check <- eau_panel_alloc %>%
  group_by(wmd_id, year, rcp, gcm) %>%
  summarise(
    total_abs  = unique(abs_density),
    sum_scaled = sum(scaled_density, na.rm = TRUE),
    diff       = sum_scaled - total_abs,
    .groups    = "drop"
  )

# Overall summary of diff
summary(alloc_check$diff)

# Look at any groups where diff is non‑negligible
alloc_check %>%
  filter(!is.na(diff), abs(diff) > 1e-6) %>%
  arrange(desc(abs(diff)))



