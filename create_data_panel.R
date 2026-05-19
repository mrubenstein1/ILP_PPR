###Construct Data Panel for Analysis ###
# Build a data panel to contain data needed for model analysis
    # takes the outputs of EAU_prop_suitable (the datatable of percent suitable 
    # habitat per year per EAU), and uses it to construct the final data table. 
    # The final data panel (i.e., final product of data import/manipulation) will be a 
    # dataframe of EAU per time period per RCP per GCM.


# ------------------------------------------------------------
# 1. Load EAU–WMD crosswalk
# ------------------------------------------------------------
eau_wmd <- read_csv("input_data/eau_wmd_lookup.csv", show_col_types = FALSE)
# eau_wmd has: eau_id, wmd_id_num, wmd_id, area_km2, x_coord, y_coord [file:2]

# ------------------------------------------------------------
# 2. Load suitable habitat tables for both RCPs
# ------------------------------------------------------------
lu_45 <- read_csv("input_data/Lu_prop_45.csv", show_col_types = FALSE)
lu_85 <- read_csv("input_data/Lu_prop_85.csv", show_col_types = FALSE)
# Columns look like rcp45_2014, rcp45_2020, ...; similarly for rcp85_* [file:1]

# ------------------------------------------------------------
# 3. Identify common years (2020–2100) across both RCP tables
# ------------------------------------------------------------
cols_45 <- grep("^rcp45_[0-9]{4}$", names(lu_45), value = TRUE)
cols_85 <- grep("^rcp85_[0-9]{4}$", names(lu_85), value = TRUE)

years_45 <- as.integer(sub("rcp45_", "", cols_45))
years_85 <- as.integer(sub("rcp85_", "", cols_85))

years_common <- sort(intersect(years_45, years_85))
years_use <- years_common[years_common >= 2020 & years_common <= 2100]

if (length(years_use) == 0) {
  stop("No overlapping years between 2020 and 2100 found in Lu_prop_45 and Lu_prop_85.")
}

cols_45_use <- paste0("rcp45_", years_use)
cols_85_use <- paste0("rcp85_", years_use)

# ------------------------------------------------------------
# 4. Build long table for RCP 4.5: EAU × year
# ------------------------------------------------------------
panel_45 <- lu_45 %>%
  select(all_of(cols_45_use)) %>%
  mutate(eau_row = row_number()) %>%
  pivot_longer(
    cols      = starts_with("rcp45_"),
    names_to  = "time_period",
    values_to = "prop_suitable"
  ) %>%
  mutate(
    year         = as.integer(str_extract(time_period, "[0-9]{4}")),
    pct_suitable = prop_suitable * 100,
    rcp          = "45"
  ) %>%
  select(eau_row, year, rcp, prop_suitable, pct_suitable)

# ------------------------------------------------------------
# 5. Build long table for RCP 8.5: EAU × year
# ------------------------------------------------------------
panel_85 <- lu_85 %>%
  select(all_of(cols_85_use)) %>%
  mutate(eau_row = row_number()) %>%
  pivot_longer(
    cols      = starts_with("rcp85_"),
    names_to  = "time_period",
    values_to = "prop_suitable"
  ) %>%
  mutate(
    year         = as.integer(str_extract(time_period, "[0-9]{4}")),
    pct_suitable = prop_suitable * 100,
    rcp          = "85"
  ) %>%
  select(eau_row, year, rcp, prop_suitable, pct_suitable)

# ------------------------------------------------------------
# 6. Stack both RCPs
# ------------------------------------------------------------
panel_suitable <- bind_rows(panel_45, panel_85)

# At this point: one row per "retained EAU row" × year × RCP [file:1]

# ------------------------------------------------------------
# 7. Attach EAU IDs and WMD IDs (crosswalk)
# ------------------------------------------------------------
# IMPORTANT: we assume that rows in lu_45 / lu_85 correspond, in order,
# to the subset of EAUs that have valid land cover data.
# we align by simple row order:

n_eau_panel <- length(unique(panel_suitable$eau_row))
if (n_eau_panel != nrow(lu_45)) {
  warning("Number of unique eau_row does not equal nrow(lu_45);",
          " please double-check alignment with eau_wmd.")
}

eau_meta <- eau_wmd %>%
  slice(seq_len(n_eau_panel)) %>%         # align by row order
  mutate(eau_row = row_number()) %>%
  select(eau_row, eau_id, wmd_id)

panel_suitable <- panel_suitable %>%
  left_join(eau_meta, by = "eau_row")

# ------------------------------------------------------------
# 8. Duplicate rows for wet vs dry GCMs
# ------------------------------------------------------------
gcm_levels <- c("wet", "dry")

panel_suitable_gcm <- panel_suitable %>%
  tidyr::crossing(gcm = gcm_levels)

# ------------------------------------------------------------
# 9. Add empty placeholders for densities, transitions, and cost
# ------------------------------------------------------------
eau_panel <- panel_suitable_gcm %>%
  mutate(
    abs_density   = NA_real_,  # to be filled later with WMD-level data
    scaled_density = NA_real_, # will depend on abs_density & prop_suitable
    trans_prob    = NA_real_,  # placeholder
    cost          = NA_real_   # placeholder
  ) %>%
  select(
    eau_id,
    wmd_id,
    year,
    rcp,
    gcm,
    prop_suitable,
    pct_suitable,
    abs_density,
    scaled_density,
    trans_prob,
    cost,
    everything()
  )

# ------------------------------------------------------------
# 10. Add 5th row per EAU-Timestep for a "stationary" scenario. 
# ------------------------------------------------------------
# Identify the key columns that define an EAU–time step
key_cols <- c("eau_id", "wmd_id", "year")

# Create a stationary row for each EAU–year
stationary_rows <- eau_panel %>%
  distinct(across(all_of(key_cols))) %>%
  mutate(
    rcp = "stationary",
    gcm = "stationary"
  )

# Make sure stationary_rows has ALL columns from eau_panel
#    For any column that's missing, add it as NA of the correct type
for (col in setdiff(names(eau_panel), names(stationary_rows))) {
  # Infer type from the existing column in eau_panel
  prototype <- eau_panel[[col]]
  if (is.numeric(prototype)) {
    stationary_rows[[col]] <- NA_real_
  } else if (is.character(prototype)) {
    stationary_rows[[col]] <- NA_character_
  } else if (is.logical(prototype)) {
    stationary_rows[[col]] <- NA
  } else {
    # Fallback: NA of same class
    stationary_rows[[col]] <- NA
  }
}

# Reorder columns to match eau_panel
stationary_rows <- stationary_rows %>%
  select(all_of(names(eau_panel)))

# Bind stationary rows to the existing panel
eau_panel <- bind_rows(
  eau_panel,
  stationary_rows
)

### Add 2020 values to percent suitable habitat for stationary at every time step. ####
    # Use the 2020 RCP 4.5 values, as this is closest to "stationary".



# 1. Extract 2020 RCP 4.5 DRY baseline per EAU
baseline_2020_45_dry <- eau_panel %>%
  filter(
    year == 2020,
    rcp  == "45",
    gcm  == "dry"
  ) %>%
  group_by(eau_id) %>%
  summarise(
    prop_suitable_baseline = first(prop_suitable),
    pct_suitable_baseline  = first(pct_suitable),
    .groups = "drop"
  )

# 2. Join baseline onto all rows, then overwrite stationary rows
eau_panel <- eau_panel %>%
  left_join(baseline_2020_45_dry, by = "eau_id") %>%
  mutate(
    prop_suitable = if_else(
      rcp == "stationary",
      prop_suitable_baseline,
      prop_suitable
    ),
    pct_suitable = if_else(
      rcp == "stationary",
      pct_suitable_baseline,
      pct_suitable
    )
  ) %>%
  select(-prop_suitable_baseline, -pct_suitable_baseline)




# ------------------------------------------------------------
# 10. Save panel for later steps
# ------------------------------------------------------------
saveRDS(eau_panel, "input_data/eau_panel.rds")
write_csv(eau_panel, "input_data/eau_panel.csv")

# ------------------------------------------------------------
# 11. Sanity check
# ------------------------------------------------------------
#### 1. Check exactly 5 rows per EAU per time period ####
rows_per_eau_year <- eau_panel %>%
  count(eau_id, year, name = "n_rows")

# How many (EAU, year) combos have each row count?
rows_per_eau_year %>%
  count(n_rows)

#### 2. Check number of time periods per EAU (should be 9) ####
years_per_eau <- eau_panel %>%
  distinct(eau_id, year) %>%
  count(eau_id, name = "n_years")

# Distribution of number of years per EAU
years_per_eau %>%
  count(n_years)


#### 3. Check total rows per EAU (should be 5 * 9 = 45) ####
rows_per_eau <- eau_panel %>%
  count(eau_id, name = "n_rows_eau")

# Distribution of total rows per EAU
rows_per_eau %>%
  count(n_rows_eau)


#### 4. Number of EAUs per WMD vs wmd_summary ####

# From the panel (using eau_id)
eaus_per_wmd_panel <- eau_panel %>%
  distinct(eau_id, wmd_id) %>%
  count(wmd_id, name = "n_eaus_panel")

eaus_per_wmd_panel

# From the original lookup (eau_wmd) – this should match wmd_summary$n_eaus
eaus_per_wmd_lookup <- eau_wmd %>%
  count(wmd_id, name = "n_eaus_lookup")

# Compare panel vs lookup
eaus_per_wmd_compare <- eaus_per_wmd_lookup %>%
  left_join(eaus_per_wmd_panel, by = "wmd_id")

eaus_per_wmd_compare

# If wmd_summary is already saved:
wmd_summary <- read_csv("input_data/wmd_summary.csv", show_col_types = FALSE)

# Compare directly to wmd_summary$n_eaus
wmd_summary_compare <- wmd_summary %>%
  select(wmd_id, n_eaus) %>%
  left_join(eaus_per_wmd_panel, by = "wmd_id")

wmd_summary_compare

#### 5. Check that for each EAU-year-RCP, the prop_suitable values are the same (no difference between wet/dry)
# For each EAU–year–RCP, count distinct prop_suitable values
suitable_check <- eau_panel %>%
  group_by(eau_id, year, rcp) %>%
  summarise(
    n_unique_prop = n_distinct(prop_suitable),
    .groups = "drop"
  )

# Summarize how many EAU–year–RCP combos have 1, 2, ... distinct values
suitable_check %>%
  count(n_unique_prop)

# Inspect any violations (where wet and dry differ)
suitable_violations <- suitable_check %>%
  filter(n_unique_prop > 1)

suitable_violations # should be 0


###. 6 For each EAU, check that stationary values are constant over time
#    (and equal to the 2020 RCP4.5 dry value)
check_stationary <- eau_panel %>%
  filter(rcp == "stationary") %>%
  group_by(eau_id) %>%
  summarise(
    n_years_stationary = n_distinct(year),
    min_pct            = min(pct_suitable, na.rm = TRUE),
    max_pct            = max(pct_suitable, na.rm = TRUE),
    .groups = "drop"
  )

check_stationary %>%
  mutate(equal_min_max = dplyr::near(min_pct, max_pct)) %>%
  count(equal_min_max)
