
##### Benefit (Duck abundance) Data ####
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

#FIRST: remove the 4WMDs we are removing from the analysis.
  #3 from Montana (Benton Lake, Bowdoin, and Northeast Montana) because the 
  #underlying hydrological data that underpinned the duck abundance estimates performed poorly 
  #in that region; and a 4th (Winsom) because of a mismatch of spatial extent with the FOREsce 
  #landcover data. 

# ------------------------------------------------------------
# Drop excluded WMDs
# ------------------------------------------------------------
wmd_exclude <- c("Windom", "Benton Lake", "Bowdoin", "Northeast Montana")

eau_panel_with_ben <- eau_panel_with_ben %>%
  filter(!wmd_id %in% wmd_exclude)

# Sanity check: confirm those WMDs are gone and row counts look right
cat("WMDs remaining:", n_distinct(eau_panel_with_ben$wmd_id), "\n")
cat("WMDs excluded:  ", 
    paste(wmd_exclude[!wmd_exclude %in% unique(eau_panel_with_ben$wmd_id)], 
          collapse = ", "), "\n")
cat("EAUs remaining:", n_distinct(eau_panel_with_ben$eau_id), "\n")


# ── Extract 2020 baseline (one row per WMD, used for 2020 and stationary) ────
ben_2020_for_join <- benefit_decadal %>%
  filter(year == 2020) %>%
  distinct(wmd_id, .keep_all = TRUE) %>%
  select(wmd_id, ben)

# ── 2020: join on wmd_id only (rcp/gcm are NA and meaningless at baseline) ───
panel_2020_with_ben <- eau_panel %>%
  filter(year == 2020) %>%
  left_join(ben_2020_for_join, by = "wmd_id") %>%
  mutate(abs_abundance = ben) %>%
  select(-ben)

# ── Future stationary: hold 2020 observed abundance constant across all years ───
panel_future_stationary_with_ben <- eau_panel %>%
  filter(year > 2020, rcp == "stationary") %>%
  left_join(ben_2020_for_join, by = "wmd_id") %>%
  mutate(abs_abundance = ben) %>%
  select(-ben)

# ── Future non-stationary: join on all four keys ──────────────────────────────
ben_future_for_join <- benefit_decadal %>%
  filter(year > 2020) %>%
  select(wmd_id, year, rcp, gcm, ben)

panel_future_nonstationary_with_ben <- eau_panel %>%
  filter(year > 2020, rcp != "stationary") %>%
  left_join(ben_future_for_join, by = c("wmd_id", "year", "rcp", "gcm")) %>%
  mutate(abs_abundance = ben) %>%
  select(-ben)

# ── Recombine ─────────────────────────────────────────────────────────────────
eau_panel_with_ben <- bind_rows(
  panel_2020_with_ben,
  panel_future_stationary_with_ben,
  panel_future_nonstationary_with_ben
) %>%
  arrange(eau_id, year, rcp, gcm)



### IV. Sanity check ####
eau_panel_with_ben %>%
  summarise(n_na = sum(is.na(abs_abundance)))

# For each WMD × year × scenario, check if all EAUs share the same abs_abundance
abs_check <- eau_panel_with_ben %>%
  group_by(wmd_id, year, rcp, gcm) %>%
  summarise(
    n_eau        = n(),
    min_abs      = min(abs_abundance, na.rm = TRUE),
    max_abs      = max(abs_abundance, na.rm = TRUE),
    n_na_abs     = sum(is.na(abs_abundance)),
    .groups = "drop"
  ) %>%
  mutate(
    all_equal = (min_abs == max_abs) | n_eau <= 1,  # or only one EAU
    any_na    = n_na_abs > 0
  )

# Look at any combinations where abs_abundance is not constant across EAUs
abs_check %>%
  filter(!all_equal & !any_na)

# ── Sanity check: abs_abundance for year 2020 matches stationary rows ───────────

# Extract 2020 baseline abundance per WMD (any rcp/gcm row will do, they're all equal)
check_2020 <- eau_panel_with_ben %>%
  filter(year == 2020) %>%
  distinct(wmd_id, abs_abundance) %>%
  rename(ben_2020 = abs_abundance)

# Extract stationary abundance per WMD (any year will do, they should all be equal)
check_stationary <- eau_panel_with_ben %>%
  filter(rcp == "stationary", gcm == "stationary") %>%
  distinct(wmd_id, abs_abundance) %>%
  rename(ben_stationary = abs_abundance)

# Compare
check_2020_vs_stationary <- check_2020 %>%
  left_join(check_stationary, by = "wmd_id") %>%
  mutate(match = near(ben_2020, ben_stationary))

n_mismatches <- sum(!check_2020_vs_stationary$match, na.rm = TRUE)
cat("WMDs where 2020 and stationary abs_abundance differ:", n_mismatches, "\n")

if (n_mismatches > 0) {
  message("⚠ Mismatched WMDs:")
  print(check_2020_vs_stationary %>% filter(!match))
} else {
  message("✓ 2020 and stationary abs_abundance values match for all WMDs.")
}




####. V. Distribute benefit by proportional habitat #####
eau_panel_alloc <- eau_panel_with_ben %>%
  group_by(wmd_id, year, rcp, gcm) %>%
  mutate(
    total_prop_suitable = sum(prop_suitable, na.rm = TRUE),
    habitat_share = ifelse(
      total_prop_suitable > 0,
      prop_suitable / total_prop_suitable,
      NA_real_
    ),
    scaled_abundance = habitat_share * abs_abundance
  ) %>%
  ungroup()

# Sanity check: scaled_density should sum back to abs_density within each group
alloc_check <- eau_panel_alloc %>%
  group_by(wmd_id, year, rcp, gcm) %>%
  summarise(
    total_abs  = first(abs_abundance),
    sum_scaled = sum(scaled_abundance, na.rm = TRUE),
    diff       = sum_scaled - total_abs,
    .groups    = "drop"
  )

summary(alloc_check$diff)

alloc_check %>%
  filter(!is.na(diff), abs(diff) > 1e-6) %>%
  arrange(desc(abs(diff)))