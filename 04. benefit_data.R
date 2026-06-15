
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

# ── Step 1: Extract the single 2020 observational anchor per WMD ─────────────
ben_2020 <- benefit_new %>%
  filter(year == 2020) %>%
  select(wmd_id, ben)

# ── Step 2: Extract future known values at 2050 and 2100 ─────────────────────
future_data <- benefit_new %>%
  filter(year %in% c(2050, 2100)) %>%
  select(wmd_id, rcp, gcm, year, ben)

# ── Step 3: Build full anchor set for interpolation ───────────────────────────
# Assign the shared 2020 observation as starting anchor for each RCP × GCM combo
anchor_2020 <- future_data %>%
  distinct(wmd_id, rcp, gcm) %>%
  left_join(ben_2020, by = "wmd_id") %>%
  mutate(year = 2020L)

anchors <- bind_rows(anchor_2020, future_data) %>%
  arrange(wmd_id, rcp, gcm, year)

# ── Step 4: Create the full decadal grid for all WMD × RCP × GCM combos ─────
decade_grid <- anchors %>%
  distinct(wmd_id, rcp, gcm) %>%
  crossing(year = as.integer(seq(2020, 2100, by = 10)))

# ── Step 5: Join anchors onto grid and linearly interpolate ───────────────────
# approx() interpolates between the three known anchor years (2020, 2050, 2100)
# treating each segment (2020→2050, 2050→2100) independently
benefit_decadal_future <- decade_grid %>%
  left_join(anchors, by = c("wmd_id", "rcp", "gcm", "year")) %>%
  arrange(wmd_id, rcp, gcm, year) %>%
  group_by(wmd_id, rcp, gcm) %>%
  mutate(
    ben = approx(
      x    = year[!is.na(ben)],
      y    = ben[!is.na(ben)],
      xout = year
    )$y
  ) %>%
  ungroup()

# ── Step 6: Bind back the single 2020 observational row ──────────────────────
ben_2020_obs <- benefit_new %>%
  filter(year == 2020) %>%
  select(wmd_id, year, rcp, gcm, ben)

benefit_decadal <- bind_rows(
  ben_2020_obs,
  benefit_decadal_future %>% filter(year > 2020)   # drop 2020 from future grid
) %>%
  arrange(wmd_id, year, rcp, gcm)

# ── Sanity check ──────────────────────────────────────────────────────────────
# Each year should have W WMDs; 2020 has 1 row/WMD, all other decades have 4
benefit_decadal %>%
  count(year, rcp, gcm)
nrow(benefit_decadal)

# 1. Inner join on all four keys
#    na_matches = "na" ensures the 2020 NA rows (rcp/gcm) join correctly
check <- benefit_new %>%
  inner_join(
    benefit_decadal,
    by         = c("wmd_id", "rcp", "gcm", "year"),
    suffix     = c("_new", "_decadal"),
    na_matches = "na"
  ) %>%
  mutate(
    diff  = ben_new - ben_decadal,
    match = near(ben_new, ben_decadal)   # near() is float-safe; avoids == rounding issues
  )

# 2. Coverage check: any benefit_new rows absent from benefit_decadal?
unmatched <- anti_join(
  benefit_new, benefit_decadal,
  by         = c("wmd_id", "rcp", "gcm", "year"),
  na_matches = "na"
)

cat("Rows in benefit_new:                        ", nrow(benefit_new),  "\n")
cat("Rows successfully matched (inner join):     ", nrow(check),        "\n")
cat("Rows in benefit_new with NO match:          ", nrow(unmatched),    "\n")

if (nrow(unmatched) > 0) {
  message("⚠ Unmatched rows:")
  print(unmatched)
}

# 3. Value check: do matched rows agree?
mismatches <- check %>% filter(!match)
cat("Value mismatches:                           ", nrow(mismatches),   "\n")

if (nrow(mismatches) > 0) {
  message("⚠ Mismatched values:")
  print(mismatches %>% select(wmd_id, year, rcp, gcm, ben_new, ben_decadal, diff))
} else {
  message("✓ All matched values are identical.")
}





#####III. Match Benefit data into EAU_Panel ####


# 1. Select only the join keys + ben from benefit_interp
benefit_for_join <- benefit_decadal %>%
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



